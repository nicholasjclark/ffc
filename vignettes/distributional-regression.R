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
library(mgcv)
library(MASS)
library(ggplot2)
library(dplyr)
theme_set(theme_bw())


## ----data-overview------------------------------------------------------------
data(mcycle, package = "MASS")
mcycle$index <- 1:NROW(mcycle)

# Check data structure
glimpse(mcycle)
cat("Total observations:", nrow(mcycle), "\n")
cat("Index range:", min(mcycle$index), "to", max(mcycle$index), "\n")


## ----acceleration-patterns, fig.cap="Head acceleration measurements showing clear heteroskedasticity over time."----
ggplot(mcycle, aes(x = index, y = accel)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_vline(
    xintercept = 110,
    linetype = "dashed",
    color = "red",
    alpha = 0.7
  ) +
  labs(
    x = "Time Index",
    y = "Head Acceleration",
    title = "Motorcycle crash data: heteroskedastic acceleration patterns"
  )


## ----variance-exploration, fig.cap="Moving variance reveals dramatic changes in acceleration variability over time."----
# Calculate moving variance with overlapping windows
window_size <- 15
mcycle$moving_var <- NA

for (i in (window_size + 1):(nrow(mcycle) - window_size)) {
  window_data <- mcycle$accel[(i - window_size):(i + window_size)]
  mcycle$moving_var[i] <- var(window_data, na.rm = TRUE)
}

ggplot(mcycle, aes(x = index)) +
  geom_line(aes(y = moving_var), color = "darkred", linewidth = 1.2) +
  geom_point(aes(y = moving_var), color = "darkred", alpha = 0.7) +
  labs(
    x = "Time Index",
    y = "Moving variance",
    title = "Time-varying variance in motorcycle crash acceleration"
  )


## ----data-split---------------------------------------------------------------
mcycle_train <- mcycle |>
  dplyr::filter(index < 110)
mcycle_test <- mcycle |>
  dplyr::filter(index >= 110)

cat("Training observations:", nrow(mcycle_train), "\n")
cat("Test observations:", nrow(mcycle_test), "\n")


## ----fit-model----------------------------------------------------------------
mcycle_gaulss_model <- ffc_gam(
  list(
    # Location parameter: time-varying mean acceleration
    accel ~ fts(index, mean_only = TRUE, time_k = 25),

    # Scale parameter: time-varying variance
    ~ fts(index, mean_only = TRUE, time_k = 25)
  ),
  family = gaulss(),
  data = mcycle_train,
  time = "index"
)


## ----model-summary------------------------------------------------------------
summary(mcycle_gaulss_model)


## ----fitted-patterns, fig.cap="Fitted location and scale parameters reveal distinct temporal patterns in crash dynamics."----
# Extract time-varying coefficients for both parameters
fts_results <- fts_coefs(mcycle_gaulss_model, summary = FALSE)

autoplot(fts_results) +
  labs(
    title = "Time-varying distributional parameters",
    x = "Time Index",
    y = "Coefficient Value"
  )


## ----extract-coefficients-----------------------------------------------------
# View coefficient structure
print(fts_results)

# Summary statistics by parameter
fts_results |>
  as_tibble() |>
  group_by(.parameter, .basis) |>
  summarise(
    mean_coef = mean(.estimate),
    sd_coef = sd(.estimate),
    .groups = "drop"
  )


## ----single-forecast----------------------------------------------------------
fc_single <- forecast(
  mcycle_gaulss_model,
  newdata = mcycle_test,
  model = "ENS"
)

# Examine forecast structure
print(fc_single)


## ----forecast-visualization, fig.cap="Distributional forecasts capture both mean trends and evolving uncertainty from the time-varying scale parameter."----
# Combine training and forecast data for visualization
forecast_plot_data <- mcycle_train |>
  as_tibble() |>
  mutate(type = "observed") |>
  bind_rows(
    mcycle_test |>
      as_tibble() |>
      bind_cols(fc_single) |>
      mutate(type = "forecast")
  )

ggplot(forecast_plot_data, aes(x = index, y = accel)) +
  geom_point(
    data = filter(forecast_plot_data, type == "observed"),
    alpha = 0.7,
    size = 2
  ) +
  geom_ribbon(
    data = filter(forecast_plot_data, type == "forecast"),
    aes(ymin = .q2.5, ymax = .q97.5),
    fill = "darkred",
    alpha = 0.2
  ) +
  geom_ribbon(
    data = filter(forecast_plot_data, type == "forecast"),
    aes(ymin = .q10, ymax = .q90),
    fill = "darkred",
    alpha = 0.3
  ) +
  geom_line(
    data = filter(forecast_plot_data, type == "forecast"),
    aes(y = .estimate),
    colour = "darkred",
    linewidth = 1.2
  ) +
  geom_point(
    data = filter(forecast_plot_data, type == "forecast"),
    color = "black",
    size = 2
  ) +
  labs(
    title = "Distributional GAM forecasts for motorcycle crash acceleration",
    x = "Time Index",
    y = "Head Acceleration"
  )


## ----coefficient-forecasting, fig.cap="Forecasted functional coefficients show parameter-specific evolution patterns with quantified uncertainty."----
# Forecast the coefficients to the test period
coef_forecast_raw <- forecast(
  fts_results,
  h = nrow(mcycle_test),
  model = "ENS"
)

# Summarize coefficient forecasts
coef_forecast <- coef_forecast_raw |>
  as_tibble() |>
  group_by(.basis, .time) |>
  summarise(
    .estimate = mean(.sim),
    .q2.5 = quantile(.sim, 0.025),
    .q97.5 = quantile(.sim, 0.975),
    .groups = "drop"
  )

# Summarize training coefficients
fts_coefs_summary <- fts_results |>
  as_tibble() |>
  group_by(.basis, .time) |>
  summarise(
    .estimate = mean(.estimate),
    .groups = "drop"
  )

# Visualize coefficient evolution and forecasts
ggplot() +
  geom_line(
    data = fts_coefs_summary,
    aes(x = .time, y = .estimate, color = .basis),
    linewidth = 1
  ) +
  geom_ribbon(
    data = coef_forecast,
    aes(x = .time, ymin = .q2.5, ymax = .q97.5, fill = .basis),
    alpha = 0.2
  ) +
  geom_line(
    data = coef_forecast,
    aes(x = .time, y = .estimate, color = .basis),
    linewidth = 1.2
  ) +
  facet_wrap(~.basis, scales = "free_y") +
  labs(
    title = "Functional coefficients: training vs forecasted",
    x = "Time Index",
    y = "Coefficient Value"
  ) +
  theme(legend.position = "none")


## ----rolling-setup------------------------------------------------------------
# Rolling forecast parameters
window_size <- 80
forecast_horizon <- 5
step_size <- 5

# Calculate forecast origins
max_start <- nrow(mcycle) - forecast_horizon
forecast_starts <- seq(1, max(1, max_start - window_size), by = step_size)

cat("Rolling forecast setup:\n")
cat("- Initial window size:", window_size, "observations\n")
cat("- Forecast horizon:", forecast_horizon, "steps\n")
cat("- Number of forecast origins:", length(forecast_starts), "\n")


## ----rolling-forecasts--------------------------------------------------------
# Storage for results
rolling_results <- list()

# Perform rolling forecasts
for (i in seq_along(forecast_starts)) {
  start_idx <- forecast_starts[i]
  end_idx <- min(start_idx + window_size - 1, nrow(mcycle))

  # Training data for this iteration
  train_data <- mcycle[start_idx:end_idx, ]

  # Test data
  test_start <- end_idx + 1
  test_end <- min(test_start + forecast_horizon - 1, nrow(mcycle))

  if (test_end <= nrow(mcycle)) {
    test_data <- mcycle[test_start:test_end, ]

    # Fit model with reduced complexity for rolling validation
    roll_model <- ffc_gam(
      list(
        accel ~ fts(
          index,
          mean_only = TRUE,
          time_k = min(15, nrow(train_data) - 5)
        ),
        ~ fts(index, mean_only = TRUE, time_k = min(15, nrow(train_data) - 5))
      ),
      family = gaulss(),
      data = train_data,
      time = "index"
    )

    # Generate forecasts
    roll_forecast <- forecast(
      roll_model,
      newdata = test_data,
      model = "ENS"
    )

    # Store results with metadata
    rolling_results[[i]] <- test_data |>
      bind_cols(roll_forecast) |>
      mutate(
        forecast_origin = i,
        train_start = start_idx,
        train_end = end_idx,
        horizon = 1:nrow(test_data)
      )
  }
}

# Combine all rolling forecast results
rolling_forecasts <- bind_rows(rolling_results)

cat("Rolling forecast evaluation completed!\n")
cat("Total forecast points:", nrow(rolling_forecasts), "\n")


## ----performance-analysis-----------------------------------------------------
# Calculate forecast errors and coverage statistics
rolling_forecasts <- rolling_forecasts |>
  mutate(
    forecast_error = accel - .estimate,
    abs_error = abs(forecast_error),
    squared_error = forecast_error^2,
    coverage_95 = (accel >= .q2.5) & (accel <= .q97.5),
    coverage_80 = (accel >= .q10) & (accel <= .q90),
    interval_width_95 = .q97.5 - .q2.5,
    interval_width_80 = .q90 - .q10
  )

# Performance by forecast horizon
horizon_stats <- rolling_forecasts |>
  group_by(horizon) |>
  summarise(
    n_forecasts = n(),
    mae = mean(abs_error, na.rm = TRUE),
    rmse = sqrt(mean(squared_error, na.rm = TRUE)),
    coverage_95 = mean(coverage_95, na.rm = TRUE),
    coverage_80 = mean(coverage_80, na.rm = TRUE),
    avg_width_95 = mean(interval_width_95, na.rm = TRUE),
    avg_width_80 = mean(interval_width_80, na.rm = TRUE),
    .groups = 'drop'
  )

print(horizon_stats)

# Overall performance summary
overall_stats <- rolling_forecasts |>
  summarise(
    total_forecasts = n(),
    overall_mae = mean(abs_error, na.rm = TRUE),
    overall_rmse = sqrt(mean(squared_error, na.rm = TRUE)),
    overall_coverage_95 = mean(coverage_95, na.rm = TRUE),
    overall_coverage_80 = mean(coverage_80, na.rm = TRUE)
  )

cat("\nOverall Rolling Forecast Performance:\n")
cat("- Total forecasts:", overall_stats$total_forecasts, "\n")
cat("- MAE:", round(overall_stats$overall_mae, 3), "\n")
cat("- RMSE:", round(overall_stats$overall_rmse, 3), "\n")
cat("- 95% coverage:", round(overall_stats$overall_coverage_95 * 100, 1), "%\n")
cat("- 80% coverage:", round(overall_stats$overall_coverage_80 * 100, 1), "%\n")


## ----rolling-visualization, fig.cap="Rolling forecast evaluation reveals consistent prediction performance across different crash phases."----
# Visualize rolling forecasts over the full time series
ggplot() +
  geom_point(data = mcycle, aes(x = index, y = accel), alpha = 0.4, size = 1) +
  geom_point(
    data = rolling_forecasts,
    aes(x = index, y = .estimate),
    color = "darkred",
    size = 1.5,
    alpha = 0.8
  ) +
  geom_errorbar(
    data = rolling_forecasts,
    aes(x = index, ymin = .q10, ymax = .q90),
    color = "darkred",
    alpha = 0.6,
    width = 0.5
  ) +
  geom_segment(
    data = rolling_forecasts,
    aes(x = index, xend = index, y = accel, yend = .estimate),
    color = "blue",
    alpha = 0.3,
    linetype = "dashed"
  ) +
  labs(
    title = "Rolling forecast evaluation across crash timeline",
    x = "Time Index",
    y = "Head Acceleration"
  )


## ----coverage-evolution, fig.cap="Cumulative coverage rates demonstrate well-calibrated prediction intervals throughout the evaluation."----
# Track coverage evolution over evaluation period
coverage_data <- rolling_forecasts |>
  arrange(index) |>
  mutate(
    running_coverage_95 = cumsum(coverage_95) / row_number(),
    running_coverage_80 = cumsum(coverage_80) / row_number()
  )

ggplot(coverage_data, aes(x = index)) +
  geom_line(aes(y = running_coverage_95), color = "darkred", linewidth = 1.2) +
  geom_line(aes(y = running_coverage_80), color = "blue", linewidth = 1.2) +
  geom_hline(
    yintercept = 0.95,
    linetype = "dashed",
    color = "darkred",
    alpha = 0.7
  ) +
  geom_hline(
    yintercept = 0.80,
    linetype = "dashed",
    color = "blue",
    alpha = 0.7
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Cumulative coverage rates validate interval calibration",
    x = "Time Index",
    y = "Cumulative Coverage Rate"
  )
