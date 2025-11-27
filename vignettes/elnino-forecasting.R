## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)


## ----setup--------------------------------------------------------------------
library(ffc)
library(ggplot2)
library(dplyr)
library(marginaleffects)
library(fable)
library(tsibble)
theme_set(theme_bw())


## ----data-overview------------------------------------------------------------
data("elnino_sst")
glimpse(elnino_sst)

# Check temporal coverage
cat("Years covered:", min(elnino_sst$year), "-", max(elnino_sst$year), "\n")
cat("Total observations:", nrow(elnino_sst), "\n")
cat("Monthly measurements per year:", elnino_sst |> 
    index_by(year) |> 
    summarise(n = n()) |> 
    pull(n) |> 
    unique(), "\n")


## ----seasonal-patterns, fig.cap="Monthly SST patterns across years."----------
ggplot(elnino_sst, aes(x = month, y = temperature, group = year, color = year)) +
  geom_line(alpha = 0.7, linewidth = 0.8) +
  scale_color_viridis_c(name = "Year") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(
    x = "Month", 
    y = "Temperature (°C)",
    title = "El Niño SST: Evolving seasonal patterns",
  )


## ----enso-identification, fig.cap="Annual mean SST anomalies highlighting major ENSO events. The 1997-1998 and 2015-2016 El Niño events stand out as extreme warm phases."----
annual_means <- elnino_sst |>
  index_by(year) |>
  summarise(
    mean_temp = mean(temperature),
    .groups = "drop"
  ) |>
  mutate(
    anomaly = mean_temp - mean(mean_temp),
    phase = case_when(
      anomaly > 0.5 ~ "El Niño",
      anomaly < -0.5 ~ "La Niña", 
      TRUE ~ "Neutral"
    )
  )

ggplot(annual_means, aes(x = year, y = anomaly, fill = phase)) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(
    values = c("El Niño" = "darkred", "La Niña" = "darkblue", "Neutral" = "gray60"),
    name = "ENSO Phase"
  ) +
  labs(
    x = "Year",
    y = "Temperature anomaly (°C)",
    title = "ENSO phases in the observational record"
  )


## ----data-split---------------------------------------------------------------
train_data <- elnino_sst |> 
  filter(year <= 2014) |>
  arrange(year, month)
test_data <- elnino_sst |> 
  filter(year > 2014) |>
  arrange(year, month)


## ----fit-model----------------------------------------------------------------
mod_elnino <- ffc_gam(
  temperature ~
    # Time-varying intercept: captures inter-annual variability
    fts(
      year, 
      mean_only = TRUE,
      time_k = 15
    ) +
    
    # Time-varying seasonal pattern: cyclic splines ensure continuity
    fts(
      month, 
      k = 12, 
      bs = "cc", 
      time_k = 15
    ),
  data = train_data,
  
  # Specify cyclic boundaries for month
  knots = list(month = c(0.5, 12.5)),
  time = "year",
  family = gaussian()
)


## ----model-summary------------------------------------------------------------
summary(mod_elnino)


## ----fitted-cycles, fig.cap="Fitted seasonal cycles for selected years. The model captures both the mean seasonal pattern and year-specific deviations."----
plot_predictions(mod_elnino, condition = c("month", "year")) +
  labs(
    x = "Month",
    y = "Temperature (°C)",
    title = "Predicted training-period seasonal cycles"
  )


## ----extract-coefficients-----------------------------------------------------
# Extract coefficients with their time series structure
func_coefs <- fts_coefs(mod_elnino, summary = FALSE, times = 10)
print(func_coefs)

# Visualize coefficient evolution
autoplot(func_coefs) +
  labs(
    title = "Time-varying basis coefficients",
  )


## ----forecast-gpdf------------------------------------------------------------
# Forecast coefficients using ETS
coef_forecast <- forecast(
  func_coefs,
  h = 4,
  times = 100,
  model = "ETS"
)

# Summarize forecasts for plotting
forecast_summary <- coef_forecast |>
  as_tibble() |>
  group_by(.basis, year) |>
  summarise(
    mean = mean(.sim),
    q05 = quantile(.sim, 0.05),
    q10 = quantile(.sim, 0.10),
    q90 = quantile(.sim, 0.90),
    q95 = quantile(.sim, 0.95),
    .groups = "drop"
  )

# Add historical values for context
historical_summary <- func_coefs |>
  as_tibble() |>
  group_by(.basis, year) |>
  summarise(
    mean = mean(.estimate),
    q05 = quantile(.estimate, 0.05),
    q10 = quantile(.estimate, 0.10),
    q90 = quantile(.estimate, 0.90),
    q95 = quantile(.estimate, 0.95),
    .groups = "drop"
  )

# Plot with ribbons for each basis function
ggplot(forecast_summary, aes(x = year)) +
  geom_ribbon(aes(ymin = q05, ymax = q95, fill = .basis), alpha = 0.2) +
  geom_ribbon(aes(ymin = q10, ymax = q90, fill = .basis), alpha = 0.3) +
  geom_line(aes(y = mean, color = .basis), linewidth = 1) +
  geom_ribbon(data = historical_summary,
              aes(ymin = q05, ymax = q95, fill = .basis), alpha = 0.2) +
  geom_ribbon(data = historical_summary,
              aes(ymin = q10, ymax = q90, fill = .basis), alpha = 0.3) +
  geom_line(data = historical_summary,
            aes(y = mean, color = .basis), linewidth = 1) +
  geom_vline(xintercept = 2014.5, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~.basis, scales = "free_y", ncol = 2) +
  scale_fill_viridis_d(guide = "none") +
  scale_color_viridis_d(guide = "none") +
  labs(
    x = "Year",
    y = "Coefficient value",
    title = "Forecasted functional basis coefficients",
  )


## ----temperature-forecast, fig.cap="Forecasted seasonal cycles for 2015-2018 with uncertainty bands."----
# Create forecast data grid
forecast_grid <- expand.grid(
  month = 1:12,
  year = 2015:2018
) |>
  as_tibble() |>
  as_tsibble(index = year, key = month) |>
  arrange(year, month)

# Generate temperature forecasts
temp_forecasts <- forecast(
  mod_elnino,
  newdata = forecast_grid,
  model = "GPDF",
  mean_model = "ETS",
  K = 3,
  type = "response",
  chains = 1
)

# Combine with actual test data for comparison
forecast_results <- temp_forecasts |>
  bind_cols(forecast_grid) |>
  left_join(
    test_data |> select(year, month, actual = temperature),
    by = c("year", "month")
  )

# Plot forecasts vs actuals
ggplot(forecast_results, aes(x = month)) +
  geom_ribbon(aes(ymin = .q2.5, ymax = .q97.5), alpha = 0.2, fill = "darkred") +
  geom_ribbon(aes(ymin = .q10, ymax = .q90), alpha = 0.3, fill = "darkred") +
  geom_line(aes(y = .estimate), color = "darkred", linewidth = 1) +
  geom_line(aes(y = actual), color = "black", linewidth = 1.2) +
  geom_point(aes(y = actual), color = "black", size = 1.5) +
  facet_wrap(~year, ncol = 2) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(
    x = "Month",
    y = "Temperature (C)",
    title = "Forecast validation: 2015-2018",
    subtitle = "Red: forecast (with 80% and 95% PIs), Black: observed"
  )


## ----spectral-validation, fig.cap="Spectral density comparison with uncertainties. The forecast preserves the key frequency components with quantified uncertainty bands."----
# Spectral analysis with uncertainty quantification
library(stats)

# Calculate spectrum for observed data
actual_ts <- ts(forecast_results$actual, frequency = 12)
actual_spec <- spectrum(actual_ts, plot = FALSE)

# Calculate spectra for forecast ensemble (using quantiles)
forecast_spectra <- list()
quantiles <- c(0.05, 0.10, 0.50, 0.90, 0.95)

for (q in quantiles) {
  forecast_q <- forecast_results[[paste0(".q", sprintf("%.0f", q * 100))]]
  if (is.null(forecast_q)) {
    # Handle different quantile naming conventions
    if (q == 0.05) forecast_q <- forecast_results$.q2.5
    else if (q == 0.95) forecast_q <- forecast_results$.q97.5
    else if (q == 0.50) forecast_q <- forecast_results$.estimate
    else next
  }
  
  forecast_ts_q <- ts(forecast_q, frequency = 12)
  forecast_spec_q <- spectrum(forecast_ts_q, plot = FALSE)
  
  forecast_spectra[[paste0("q", sprintf("%.0f", q * 100))]] <- data.frame(
    frequency = forecast_spec_q$freq,
    power = forecast_spec_q$spec,
    quantile = q
  )
}

# Combine spectral estimates
forecast_spec_df <- do.call(rbind, forecast_spectra)

# Calculate summary statistics across quantiles for each frequency
forecast_summary <- forecast_spec_df |>
  group_by(frequency) |>
  summarise(
    median = median(power),
    q05 = quantile(power, 0.05, na.rm = TRUE),
    q10 = quantile(power, 0.10, na.rm = TRUE), 
    q90 = quantile(power, 0.90, na.rm = TRUE),
    q95 = quantile(power, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

# Observed spectrum data
observed_df <- data.frame(
  frequency = actual_spec$freq,
  power = actual_spec$spec
)

# Create plot with uncertainty bands
ggplot() +
  # Forecast uncertainty bands
  geom_ribbon(
    data = forecast_summary,
    aes(x = frequency, ymin = log(q05), ymax = log(q95)),
    alpha = 0.2, fill = "darkred"
  ) +
  geom_ribbon(
    data = forecast_summary, 
    aes(x = frequency, ymin = log(q10), ymax = log(q90)),
    alpha = 0.3, fill = "darkred"
  ) +
  # Forecast median
  geom_line(
    data = forecast_summary,
    aes(x = frequency, y = log(median)),
    color = "darkred", linewidth = 1.2
  ) +
  # Observed spectrum
  geom_line(
    data = observed_df,
    aes(x = frequency, y = log(power)),
    color = "black", linewidth = 1.2
  ) +
  # Highlight annual frequency
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.6) +
  annotate("text", x = 1.1, y = max(log(observed_df$power)) * 0.9, 
           label = "Annual cycle", hjust = 0, size = 3) +
  labs(
    x = "Frequency (cycles per year)",
    y = "Log spectral density",
    title = "Spectral validation with forecast uncertainty",
    subtitle = "Black: observed, Red: forecast median with 80% and 90% bands"
  ) +
  theme(legend.position = "none")


## ----fable-integration, fig.height=9, fig.cap="El Niño SST forecasts using ensemble model with fable integration. Red ribbons show prediction intervals, black lines show observed values."----
# Convert ffc forecasts to fable format
fable_forecasts <- as_fable(
  mod_elnino,
  newdata = test_data,
  model = "ENS",
  mean_model = "ETS",
  type = "response",
  chains = 1
)

# Evaluate forecast accuracy using fable metrics
accuracy(fable_forecasts, test_data)

# Create publication-ready plots using fable's autoplot
fable_forecasts |>
  autoplot(data = train_data, show_gap = FALSE, fill = "darkred") +
  labs(
    x = "Year",
    y = "Temperature (C)",
    title = "Forecasts broken down by month"
  )

