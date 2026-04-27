library(tidyverse)
library(duckdbfs)
library(duckdb)

# set variables 
my_model_id <- "LogReg_Bloom"
forecast_horizon <- 30
n_members <- 300

# reading in target data
url <- "https://amnh1.osn.mghpcc.org/bio230121-bucket01/vera4cast/targets/project_id=vera4cast/duration=P1D/daily-insitu-targets.csv.gz"
targets <- read_csv(url, show_col_types = FALSE)
targets <- targets %>% filter(variable == "Bloom_binary_mean")

# past weather data
weather_stage3 <- vera4castHelpers::noaa_stage3()
df_historical <- weather_stage3 |> 
  dplyr::filter(datetime >= ymd('2022-07-06'),
                site_id %in% unique(targets$site_id))|>
  dplyr::collect()
# processing past weather data 
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
# setting up dates 
forecast_date <- Sys.Date() 
noaa_date <- forecast_date - days(1)
forecasted_dates <- seq(from = ymd(forecast_date), to = ymd(forecast_date) + forecast_horizon, by = "day")


# driver weather data 
weather <- vera4castHelpers::noaa_stage2(start_date = as.character(noaa_date))
df_future <- weather |> 
  dplyr::filter(site_id %in% unique(targets$site_id)) %>%
  dplyr::collect()

weather_future_daily <- df_future |> 
  mutate(datetime = as_date(datetime)) |> 
  # mean daily forecasts at each site per ensemble
  group_by(datetime, site_id, parameter, variable) |> 
  summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
  # convert air temperature to Celsius if it is included in the weather data
  mutate(prediction = ifelse(variable == "air_temperature", prediction - 273.15, prediction)) |> 
  pivot_wider(names_from = variable, values_from = prediction) 
weather_future_daily <- weather_future_daily%>%
  mutate(cumul_wind = abs(northward_wind) + abs(eastward_wind),
         high_rad = case_when(
           surface_downwelling_shortwave_flux_in_air < 150 ~ 0,
           surface_downwelling_shortwave_flux_in_air >=150 ~ 1
         ))
# remodeling our target data 
targets_model <- targets |> 
  pivot_wider(names_from = 'variable', values_from = 'observation') |> 
  left_join(weather_past_daily, 
            by = c("datetime","site_id"))
targets_model$Bloom_binary_mean <- as.factor(targets_model$Bloom_binary_mean)
targets_model$high_rad <- as.factor(targets_model$high_rad)


##########################################################################################
# starting the loops 
forecast_df <- NULL
ens_n <- 1

# loop through both sites 
for(i in 1:2) {  
     #i = 1 
  #print(i)}
  curr_site <- unique(targets_model$site_id)[i]
  
  site_target <- targets_model |>
    filter(site_id == curr_site) %>%
    na.omit()
  
  noaa_future_site <- weather_future_daily |> 
    filter(site_id == curr_site)
  
  weather_ensemble_names <- unique(noaa_future_site$parameter)
  
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
  
  #looping through dates 
  for (t in 1:length(forecasted_dates)) {
    #print(t)
    # looping through infinite ensemble members 
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
                             variable = "Bloom_binary_mean") 
       
      forecast_df <- dplyr::bind_rows(forecast_df, curr_site_df)
      
    }
  }
  message(curr_site, ' forecast run')
}

# reformatting our outputs
forecast_long <- forecast_df %>% group_by(site_id, datetime) %>%
  summarize(mean = mean(prediction, na.rm=TRUE),
            sd = sd(prediction))%>%
  pivot_longer(cols = c("mean", "sd" ),
               names_to = "parameter",
               values_to = "prediction")

forecast_df_EFI <- forecast_long %>%
  filter(datetime > forecast_date) %>%
  mutate(model_id = my_model_id,
         reference_datetime = forecast_date,
         family = 'bernoulli',
         duration = 'P1D',
         depth_m = case_when(site_id == "bvre" ~ 1.5,
                   site_id == "fcre" ~ 1.6
           
         ),
         variable = "Bloom_binary_mean",
         project_id = 'vera4cast') %>%
  select(datetime, reference_datetime, duration, site_id, family, depth_m, parameter, variable, prediction, model_id, project_id)

#---------------------------#

theme <- 'biological'
date <- forecast_df_EFI$reference_datetime[1]
forecast_file <- paste(theme, date, forecast_name, sep = '-')
forecast_name <- paste0(forecast_df_EFI$model_id[1], ".csv")
forecast_file <- paste(theme, date, forecast_name, sep = '-')

write_csv(forecast_df_EFI, forecast_file)

vera4castHelpers::forecast_output_validator(forecast_file)
vera4castHelpers:::submit(forecast_file)
