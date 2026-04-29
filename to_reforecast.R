library(tidyverse)
library(lubridate)
library(zoo)
library(duckdbfs)
library(duckdb)
library(scoringRules)

my_model_id <- "LogReg_Bloom"
forecast_horizon <- 30
n_members <- 300

url <- "https://amnh1.osn.mghpcc.org/bio230121-bucket01/vera4cast/targets/project_id=vera4cast/duration=P1D/daily-insitu-targets.csv.gz"

targets <- read_csv(url, show_col_types = FALSE) %>%
  filter(variable == "Bloom_binary_mean") 

weather_stage3 <- vera4castHelpers::noaa_stage3()

df_historical <- weather_stage3 |>
  filter(site_id %in% unique(targets$site_id),
         datetime < as.Date("2025-01-01")) |>
  collect()

weather_past_daily <- df_historical  |> 
  mutate(datetime = as_date(datetime)) |> 
  group_by(datetime, site_id, variable) |> 
  summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
  # convert air temperature to Celsius if it is included in the weather data
  mutate(prediction = ifelse(variable == "air_temperature", prediction - 273.15, prediction)) |> 
  pivot_wider(names_from = variable, values_from = prediction)
weather_past_daily <- weather_past_daily %>%
  mutate(cumul_wind = abs(northward_wind) + abs(eastward_wind),
         high_rad = case_when(
           surface_downwelling_shortwave_flux_in_air < 150 ~ 0,
           surface_downwelling_shortwave_flux_in_air >=150 ~ 1
         ))
targets_model <- targets |> 
  filter(datetime < as.Date("2025-01-01"))|>
  pivot_wider(names_from = 'variable', values_from = 'observation') |> 
  left_join(weather_past_daily, 
            by = c("datetime","site_id"))
targets_model$Bloom_binary_mean <- as.factor(targets_model$Bloom_binary_mean)
targets_model$high_rad <- as.factor(targets_model$high_rad)



dates <- seq(
  as.Date("2025-01-01"),
  as.Date("2026-01-01"),
  by = "7 days"
)

forecast_df <- NULL
ens_n <- 1
for(i in 1:2) {  
  #i = 1 
  #print(i)}
  curr_site <- unique(targets_model$site_id)[i]
  print(curr_site)
  site_target <- targets_model |>
    filter(site_id == curr_site) %>%
    na.omit()
  
  fit <- glm(Bloom_binary_mean ~ cumul_wind + surface_downwelling_longwave_flux_in_air+factor(high_rad)
             ,
             data = site_target,
             family = "binomial"
  )
  coeffs <- round(fit$coefficients, 2)
  #coeffs
  fit_summary <- summary(fit)
  #sigma for process uncertainty
  #### NEED 
  sigma <- sd(resid(fit, type = "response"), na.rm=TRUE)
  
  # parameters for parameter uncertainty 
  params_se <- fit_summary$coefficients[,2]
  param_df <- data.frame(beta1 = rnorm(n_members, coeffs[1], params_se[1]),
                         beta2 = rnorm(n_members, coeffs[2], params_se[2]),
                         beta3 = rnorm(n_members, coeffs[3], params_se[3]),
                         beta4 = rnorm(n_members, coeffs[4], params_se[4]))
  
  for (d in dates) {
    
  forecast_date <- as.Date(d)
  print(forecast_date)
  noaa_date <- forecast_date - days(1)
  
  weather_future_s3 <- vera4castHelpers::noaa_stage2(
    start_date = as.character(noaa_date)
  )
  
  df_future <- weather_future_s3 |> 
    filter(site_id == curr_site) %>%
    collect()
  
  weather_future_daily <- df_future |> 
    mutate(datetime = as_date(datetime)) |> 
    group_by(datetime, site_id, parameter, variable) |> 
    summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
    mutate(prediction = ifelse(variable == "air_temperature",
                               prediction - 273.15,
                               prediction)) |> 
    pivot_wider(names_from = variable, values_from = prediction)
  weather_future_daily <- weather_future_daily%>%
    mutate(cumul_wind = abs(northward_wind) + abs(eastward_wind),
           high_rad = case_when(
             surface_downwelling_shortwave_flux_in_air < 150 ~ 0,
             surface_downwelling_shortwave_flux_in_air >=150 ~ 1
           ))
  
  noaa_future_site <- weather_future_daily |> 
    filter(site_id == curr_site)
  
  weather_ensemble_names <- unique(noaa_future_site$parameter)
  
  forecasted_dates <- seq(
    from = ymd(forecast_date),
    to = ymd(forecast_date) + forecast_horizon,
    by = "day"
  )
  for (t in 1:length(forecasted_dates)) {
    
    for(ens in 1:n_members){
      if (ens_n <= 29) {
        ens_n <- ens_n+1}
      else {ens_n = 0}
      
      met_ens <- weather_ensemble_names[ens_n]
      
      #pull driver ensemble for the relevant date
      temp_driv <- weather_future_daily %>%
        filter(datetime == forecasted_dates[t],
               site_id == curr_site,
               parameter == ens_n)
      if (nrow(temp_driv) == 0){
        next
      }
      # tweaking model parameters 
      fit$coefficients[1] <- param_df$beta1[ens]
      fit$coefficients[2] <- param_df$beta2[ens]
      fit$coefficients[3] <- param_df$beta3[ens]
      fit$coefficients[4] <- param_df$beta4[ens]
      probabilities <- predict(fit, newdata = temp_driv, type = "response")
      
      
      threshold <- 1
      # adding process uncertainty but keeping our probabilities within normal range 
      new_probs <- pmin(threshold, pmax(0, probabilities[1] + rnorm(n=1, 0, sigma)))
      
      curr_site_df <- tibble(datetime = forecasted_dates[t],
                             site_id = curr_site,
                             parameter = ens,
                             prediction = new_probs,
                             reference_datetime = forecast_date,
                             variable = "Bloom_binary_mean") 
      
      forecast_df <- dplyr::bind_rows(forecast_df, curr_site_df)
      
    }
  }
  }
}
write.csv(forecast_df, "reforecast.csv")

#############################################
# Plotting Section
#############################################

new_targets <- targets |> 
  filter(datetime >= as.Date("2025-01-01"))|>
  pivot_wider(names_from = 'variable', values_from = 'observation') 
write.csv(new_targets, "new_targets.csv")

summary_df <- forecast_df %>% group_by(datetime, reference_datetime, site_id) %>%
  summarize(mean_pred = mean(prediction),
            sd_pred = sd(prediction))%>%
  left_join(new_targets, by = c("datetime", "site_id"))%>%
  mutate(crps = crps_norm(Bloom_binary_mean, mean_pred, sd_pred),
         site_id = factor(site_id, levels = c("bvre", "fcre"), labels = c("Beaverdam", "Falling Creek")))%>%
  group_by(datetime, site_id)%>%
  summarize(mean_crps = mean(crps))


all_results <- arrow::open_dataset("s3://anonymous@bio230121-bucket01/vera4cast/forecasts/bundled-summaries/project_id=vera4cast/duration=P1D/variable=Bloom_binary_mean/model_id=climatology?endpoint_override=amnh1.osn.mghpcc.org")
climatology <- all_results |> filter(datetime >= as.Date("2025-01-01"),
                            datetime <= as.Date("2026-01-01"))%>%
  dplyr::collect()

clim_comp <- climatology %>% left_join(new_targets, by = c("datetime", "site_id") ) %>%
  mutate(crps = crps_norm(Bloom_binary_mean, mean = mean, sd = sd),
         site_id = factor(site_id, levels = c("bvre", "fcre"), labels = c("Beaverdam", "Falling Creek")))
  


all_results <- arrow::open_dataset("s3://anonymous@bio230121-bucket01/vera4cast/forecasts/bundled-summaries/project_id=vera4cast/duration=P1D/variable=Bloom_binary_mean/model_id=persistenceRW?endpoint_override=amnh1.osn.mghpcc.org")
persistence <- all_results |> filter(datetime >= as.Date("2025-01-01"),
                                     datetime <= as.Date("2026-01-01"))%>%
  dplyr::collect()

pers_comp <- persistence %>% left_join(new_targets, by = c("datetime", "site_id") ) %>%
  mutate(crps = crps_norm(Bloom_binary_mean, mean = mean, sd = sd),
         site_id = factor(site_id, levels = c("bvre", "fcre"), labels = c("Beaverdam", "Falling Creek")))

clim_datetime <- clim_comp %>% group_by(datetime, site_id, model_id)%>%
  summarize(mean_crps = mean(crps))
pers_datetime <- pers_comp %>% group_by(datetime, site_id, model_id)%>%
  summarize(mean_crps = mean(crps))

ggplot() + 
  geom_line(data = pers_datetime, aes(x=as.Date(datetime), y=mean_crps, color = "Persistence")) +
  geom_line(data = clim_datetime, aes(x=as.Date(datetime), y = mean_crps, color = "Climatology"))+
  geom_line(data = summary_df, aes(x=as.Date(datetime), y=mean_crps, color = "My Model"))+
  scale_color_manual(values = c( "slateblue", "#006633", "red"))+
  facet_wrap(~site_id, ncol = 1)+labs(x="Date", y = "Mean CRPS", title = "CRPS for Three Bloom Forecasts in 2025")+ 
  theme(legend.title = element_blank())

summary_df <- forecast_df %>% group_by(datetime, reference_datetime, site_id) %>%
  summarize(mean_pred = mean(prediction),
            sd_pred = sd(prediction))%>%
  left_join(new_targets, by = c("datetime", "site_id"))%>%
  na.omit(cols = "Bloom_binary_mean")%>%
  mutate(crps = crps_norm(Bloom_binary_mean, mean_pred, sd_pred),
         site_id = factor(site_id, levels = c("bvre", "fcre"), labels = c("Beaverdam", "Falling Creek")),
         horizon = as.integer(as.Date(datetime)-as.Date(reference_datetime)))%>%
  group_by(horizon, site_id)%>%
  summarize(mean_crps = mean(crps))



clim_horizon <- clim_comp %>% mutate(horizon = as.integer(datetime-reference_datetime))%>%
  na.omit(cols = "Bloom_binary_mean")%>%
  group_by(horizon, site_id) %>%
  summarize(mean_crps = mean(crps))
pers_horizon <- pers_comp  %>% mutate(horizon = as.integer(datetime-reference_datetime))%>%
  na.omit(cols = "Bloom_binary_mean")%>%
  group_by(horizon, site_id) %>%
  summarize(mean_crps = mean(crps))

ggplot() + 
  geom_line(data = pers_horizon, aes(x=horizon, y=mean_crps, color = "Persistence")) +
  geom_line(data = clim_horizon, aes(x=horizon, y = mean_crps, color = "Climatology"))+
  geom_line(data = summary_df, aes(x=horizon, y=mean_crps, color = "My Model"))+
  scale_color_manual(values = c("slateblue", "#006633", "red"))+
  facet_wrap(~site_id, ncol = 1)+labs(x="Horizon (Days)", y = "Mean CRPS", title = "CRPS for Three Bloom Forecasts across Horizons")+ 
  theme(legend.title = element_blank())
