# Distributional Regression Vignette Script
# Demonstrates motorcycle crash data with time-varying location and scale using gaulss()
# Using rolling forecast evaluation to assess model performance

# Load packages
devtools::load_all()
library(mgcv)
library(MASS)
library(ggplot2)
library(dplyr)

# Set up plotting theme
theme_set(theme_bw())

# =============================================================================
# DATA PREPARATION: Motorcycle Crash Data
# =============================================================================

# Load and prepare motorcycle crash data
data(mcycle, package = "MASS")
mcycle$index <- 1:NROW(mcycle)

cat("=== Motorcycle Crash Data Overview ===\n")
cat("Total observations:", nrow(mcycle), "\n")
cat("Index range:", min(mcycle$index), "to", max(mcycle$index), "\n")
cat("Acceleration summary:\n")
print(summary(mcycle$accel))

# Initial split for model development
mcycle_train <- mcycle |>
  dplyr::filter(index < 110)
mcycle_test <- mcycle |>
  dplyr::filter(index >= 110)

# Visualize the data
p1 <- ggplot(mcycle, aes(x = index, y = accel)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_vline(xintercept = 110, linetype = "dashed", color = "red", alpha = 0.7) +
  labs(title = "Motorcycle Crash Data: Head Acceleration vs Time",
       subtitle = "Red line shows initial train/test split",
       x = "Time Index", y = "Head Acceleration") +
  theme_minimal()

print(p1)

# =============================================================================
# DISTRIBUTIONAL REGRESSION MODEL FITTING
# =============================================================================

cat("\n=== Fitting Distributional GAM Model ===\n")
cat("Using gaulss() family for Gaussian location-scale modeling\n")

# Fit distributional regression model with time-varying location and scale
mcycle_gaulss_model <- ffc_gam(
  list(
    accel ~ fts(index, mean_only = TRUE, time_k = 25),
    ~ fts(index, mean_only = TRUE, time_k = 25)
  ),
  family = gaulss(),
  data = mcycle_train,
  time = "index"
)

# Model summary
cat("\nModel fitted successfully!\n")
summary(mcycle_gaulss_model)

# =============================================================================
# MODEL INSPECTION AND DIAGNOSTICS
# =============================================================================

cat("\n=== Model Inspection ===\n")

# Extract and plot time-varying coefficients for both parameters
gratia::draw(mcycle_gaulss_model)

# Plot functional time series coefficients
fts_coefs(mcycle_gaulss_model, summary = FALSE) |> autoplot()

# =============================================================================
# SINGLE FORECAST DEMONSTRATION
# =============================================================================

cat("\n=== Single Forecast Demonstration ===\n")

# Generate forecasts for test period
cat("Generating forecasts for", nrow(mcycle_test), "time points ahead...\n")
fc_single <- forecast(
  mcycle_gaulss_model,
  newdata = mcycle_test,
  model = "ENS"
)

# Visualize single forecast
forecast_plot_data <- mcycle_train |>
  as_tibble() |>
  mutate(type = "observed") |>
  bind_rows(
    mcycle_test |>
      as_tibble() |>
      bind_cols(fc_single) |>
      mutate(type = "forecast")
  )

p2 <- ggplot(forecast_plot_data, aes(x = index, y = accel)) +
  geom_point(data = filter(forecast_plot_data, type == "observed"),
             alpha = 0.7, size = 2) +
  geom_ribbon(data = filter(forecast_plot_data, type == "forecast"),
              aes(ymin = .q2.5, ymax = .q97.5),
              fill = "darkred", alpha = 0.2) +
  geom_ribbon(data = filter(forecast_plot_data, type == "forecast"),
              aes(ymin = .q10, ymax = .q90),
              fill = "darkred", alpha = 0.3) +
  geom_line(data = filter(forecast_plot_data, type == "forecast"),
            aes(y = .estimate),
            colour = "darkred", size = 1.2) +
  geom_point(data = filter(forecast_plot_data, type == "forecast"),
             color = "black", size = 2) +
  labs(title = "Motorcycle Crash Data: Distributional GAM Forecast",
       subtitle = "Training data (points) and out-of-sample forecasts (red)",
       x = "Time Index", y = "Head Acceleration") +
  theme_minimal()

print(p2)

# =============================================================================
# COEFFICIENT FORECAST INSPECTION
# =============================================================================

# Extract and forecast the functional coefficients
fts_coefs_training <- fts_coefs(mcycle_gaulss_model, summary = FALSE)

# Forecast the coefficients themselves to the test period
coef_forecast_raw <- forecast(
  fts_coefs_training,
  h = nrow(mcycle_test),
  model = "ENS"
)

# Summarize forecasts
coef_forecast <- coef_forecast_raw |>
  as_tibble() |>
  group_by(.basis, .time) |>
  summarise(
    .estimate = mean(.sim),
    .q2.5 = quantile(.sim, 0.025),
    .q97.5 = quantile(.sim, 0.975),
    .groups = "drop"
  )

# Summarize training coefficients for consistent plotting
fts_coefs_summary <- fts_coefs_training |>
  as_tibble() |>
  group_by(.basis, .time) |>
  summarise(
    .estimate = mean(.estimate),
    .groups = "drop"
  )

# Visualize training coefficients vs forecasted coefficients
p_coefs <- ggplot() +
  # Training coefficients
  geom_line(data = fts_coefs_summary,
            aes(x = .time, y = .estimate, color = .basis),
            size = 1) +
  # Forecasted coefficients
  geom_ribbon(data = coef_forecast,
              aes(x = .time, ymin = .q2.5, ymax = .q97.5, fill = .basis),
              alpha = 0.2) +
  geom_line(data = coef_forecast,
            aes(x = .time, y = .estimate, color = .basis),
            size = 1.2) +
  facet_wrap(~ .basis, scales = "free_y") +
  labs(title = "Functional Coefficients: Training vs Forecasted",
       x = "Time Index", y = "Coefficient Value") +
  theme_minimal() +
  theme(legend.position = "none")

print(p_coefs)

# =============================================================================
# ROLLING FORECAST EVALUATION
# =============================================================================

cat("\n=== Rolling Forecast Evaluation ===\n")

# Set up rolling forecast parameters
window_size <- 80    # Initial training window
min_window <- 80     # Minimum training window
forecast_horizon <- 5 # Forecast 5 steps ahead each time
step_size <- 5       # Move window by 5 observations each time

# Calculate the rolling forecast range
max_start <- nrow(mcycle) - forecast_horizon
forecast_starts <- seq(1, max(1, max_start - window_size), by = step_size)

cat("Rolling forecast setup:\n")
cat("- Initial window size:", window_size, "observations\n")
cat("- Forecast horizon:", forecast_horizon, "steps\n")
cat("- Step size:", step_size, "observations\n")
cat("- Number of forecast origins:", length(forecast_starts), "\n")

# Storage for rolling forecasts
rolling_results <- list()

# Perform rolling forecasts
for (i in seq_along(forecast_starts)) {
  start_idx <- forecast_starts[i]
  end_idx <- min(start_idx + window_size - 1, nrow(mcycle))

  # Training data for this iteration
  train_data <- mcycle[start_idx:end_idx, ]

  # Test data (next forecast_horizon observations)
  test_start <- end_idx + 1
  test_end <- min(test_start + forecast_horizon - 1, nrow(mcycle))

  if (test_end <= nrow(mcycle)) {
    test_data <- mcycle[test_start:test_end, ]

    # Fit model on training data
    roll_model <- ffc_gam(
      list(
        accel ~ fts(index, mean_only = TRUE, time_k = min(15, nrow(train_data) - 5)),
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

    # Store results
    rolling_results[[i]] <- test_data |>
      bind_cols(roll_forecast) |>
      mutate(
        forecast_origin = i,
        train_start = start_idx,
        train_end = end_idx,
        horizon = 1:nrow(test_data)
      )

    cat("Completed forecast", i, "of", length(forecast_starts),
        "- training:", start_idx, "to", end_idx,
        "- testing:", test_start, "to", test_end, "\n")
  }
}

# Combine rolling forecast results
rolling_forecasts <- bind_rows(rolling_results)

cat("\nRolling forecast evaluation completed!\n")
cat("Total forecast points:", nrow(rolling_forecasts), "\n")

# =============================================================================
# ROLLING FORECAST ANALYSIS
# =============================================================================

cat("\n=== Rolling Forecast Performance Analysis ===\n")

# Calculate forecast errors
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

# Summary statistics by forecast horizon
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

print("Forecast Performance by Horizon:")
print(horizon_stats)

# Overall performance
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

# =============================================================================
# ROLLING FORECAST VISUALIZATION
# =============================================================================

# Plot rolling forecasts over time
p3 <- ggplot() +
  # Original data points
  geom_point(data = mcycle, aes(x = index, y = accel),
             alpha = 0.4, size = 1) +
  # Rolling forecast points and intervals
  geom_point(data = rolling_forecasts,
             aes(x = index, y = .estimate),
             color = "darkred", size = 1.5, alpha = 0.8) +
  geom_errorbar(data = rolling_forecasts,
                aes(x = index, ymin = .q10, ymax = .q90),
                color = "darkred", alpha = 0.6, width = 0.5) +
  geom_segment(data = rolling_forecasts,
               aes(x = index, xend = index, y = accel, yend = .estimate),
               color = "blue", alpha = 0.3, linetype = "dashed") +
  labs(title = "Rolling Forecast Evaluation: Motorcycle Crash Data",
       subtitle = paste("Red points: forecasts | Blue dashes: forecast errors |",
                        "Error bars: 80% prediction intervals"),
       x = "Time Index", y = "Head Acceleration") +
  theme_minimal()

print(p3)

# Forecast error evolution
p4 <- ggplot(rolling_forecasts, aes(x = index, y = forecast_error)) +
  geom_point(aes(color = factor(horizon)), size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_color_viridis_d(name = "Forecast\nHorizon") +
  labs(title = "Forecast Errors Over Time",
       subtitle = "Rolling forecast evaluation results",
       x = "Time Index", y = "Forecast Error (Actual - Predicted)") +
  theme_minimal()

print(p4)

# Coverage evolution
coverage_data <- rolling_forecasts |>
  arrange(index) |>
  mutate(
    running_coverage_95 = cumsum(coverage_95) / row_number(),
    running_coverage_80 = cumsum(coverage_80) / row_number()
  )

p5 <- ggplot(coverage_data, aes(x = index)) +
  geom_line(aes(y = running_coverage_95), color = "red", size = 1.2) +
  geom_line(aes(y = running_coverage_80), color = "blue", size = 1.2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "blue", alpha = 0.7) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Cumulative Coverage Rates",
       subtitle = "Red: 95% intervals | Blue: 80% intervals | Dashed: nominal rates",
       x = "Time Index", y = "Cumulative Coverage Rate") +
  theme_minimal()

print(p5)
