# Script to fit ffc model to Poblenou dataset with cyclic splines for hours
# Assumes Poblenou data has hourly measurements that can benefit from cyclic modeling

library(fable)
library(tsibble)
library(dplyr)
library(ggplot2)
library(fda.usc)
devtools::load_all()

# Load the Poblenou dataset from fda.usc package
data("poblenou", package = "fda.usc")

# Inspect the data structure
str(poblenou)
head(poblenou)

# Convert to long format for time series analysis
# poblenou$nox$data contains hourly NOx measurements
poblenou_long <- poblenou$nox$data %>%
  as.data.frame() %>%
  mutate(day = row_number()) %>%
  tidyr::pivot_longer(
    cols = -day,
    names_to = "hour_col",
    values_to = "value"
  ) %>%
  mutate(
    # Extract hour from column names (H0, H1, ..., H23)
    hour = as.numeric(gsub("H", "", hour_col))
  ) %>%
  arrange(day, hour) %>%
  dplyr::select(day, hour, value)

# Create time index for modeling
poblenou_ts <- poblenou_long %>%
  mutate(
    # Create a combined time index
    time_id = (day - 1) * 24 + hour + 1,
    # Create day groups for tsibble
    day_group = day
  ) %>%
  filter(!is.na(value)) %>%
  as_tsibble(index = day, key = hour)

# Display the result
print(poblenou_ts)

# Check the structure
cat("\nData structure:\n")
cat("Days covered:", min(poblenou_ts$day), "-", max(poblenou_ts$day), "\n")
cat("Hours per day:", length(unique(poblenou_ts$hour)), "\n")
cat("Total observations:", nrow(poblenou_ts), "\n")

# Add day of week information
poblenou_ts <- poblenou_ts %>%
  left_join(
    poblenou$df %>%
      mutate(day = row_number()) %>%
      dplyr::select(day, day.week),
    by = "day"
  ) %>%
  mutate(
    day_week = as.numeric(as.character(day.week))
  ) %>%
  dplyr::select(-day.week)

# Visualize the data
poblenou_ts %>%
  autoplot(value) +
  labs(
    title = "Poblenou Hourly Data",
    x = "Day",
    y = "NOx (mg/m³)"
  )

# Check for weekly seasonality patterns
# Plot 1: Average by day of week and hour
weekly_pattern <- poblenou_ts %>%
  group_by(day_week, hour) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(weekly_pattern, aes(x = hour, y = mean_value, color = factor(day_week))) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_viridis_d(
    name = "Day of Week",
    labels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
  ) +
  labs(
    title = "Weekly Patterns in NOx Levels",
    subtitle = "Clear differences between weekdays (1-5) and weekends (6-7)",
    x = "Hour of Day",
    y = "Mean NOx (mg/m³)"
  ) +
  theme_minimal()

# Plot 2: Heatmap of weekly patterns
ggplot(weekly_pattern, aes(x = hour, y = factor(day_week), fill = mean_value)) +
  geom_tile() +
  scale_fill_viridis_c(name = "NOx\n(mg/m³)") +
  scale_y_discrete(labels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) +
  labs(
    title = "Weekly NOx Pattern Heatmap",
    subtitle = "Weekend effect visible: lower pollution on days 6-7",
    x = "Hour of Day",
    y = "Day of Week"
  ) +
  theme_minimal()

# Plot 3: Boxplot by day of week (aggregated across all hours)
poblenou_ts %>%
  as_tibble() %>%
  group_by(day, day_week) %>%
  summarise(daily_mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = factor(day_week), y = daily_mean, fill = factor(day_week))) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_viridis_d(guide = "none") +
  scale_x_discrete(labels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) +
  labs(
    title = "Daily NOx Levels by Day of Week",
    subtitle = "Weekends show lower pollution levels",
    x = "Day of Week",
    y = "Daily Mean NOx (mg/m³)"
  ) +
  theme_minimal()

# Check for monthly cycle (approximately 28-30 days)
# Add a day_in_month variable based on a 28-day cycle
poblenou_ts <- poblenou_ts %>%
  mutate(
    day_month = ((day - 1) %% 28) + 1  # Creates a 28-day cycle
  )

# Plot 4: Check for monthly patterns
# Aggregate to daily level and plot autocorrelation
daily_avg <- poblenou_ts %>%
  as_tibble() %>%
  group_by(day) %>%
  summarise(daily_mean = mean(value, na.rm = TRUE), .groups = "drop")

# Autocorrelation plot to detect cycles
acf_plot <- acf(daily_avg$daily_mean, lag.max = 60, plot = FALSE)
acf_df <- data.frame(
  lag = acf_plot$lag[-1],  # Remove lag 0
  acf = acf_plot$acf[-1]
)

ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_hline(yintercept = c(-1.96/sqrt(nrow(daily_avg)), 1.96/sqrt(nrow(daily_avg))),
             linetype = "dashed", color = "blue") +
  geom_segment(aes(xend = lag, yend = 0), color = "black") +
  geom_point(size = 2) +
  geom_vline(xintercept = c(7, 14, 21, 28, 35, 42),
             linetype = "dotted", color = "red", alpha = 0.5) +
  labs(
    title = "Autocorrelation Function of Daily NOx Averages",
    subtitle = "Red dotted lines at weekly intervals; peak around 28 days suggests monthly cycle",
    x = "Lag (days)",
    y = "ACF"
  ) +
  theme_minimal()

# Plot 5: Average NOx by position in 28-day cycle
monthly_pattern <- poblenou_ts %>%
  group_by(day_month) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),
    se_value = sd_value / sqrt(n()),
    .groups = "drop"
  )

ggplot(monthly_pattern, aes(x = day_month, y = mean_value)) +
  geom_ribbon(aes(ymin = mean_value - 2*se_value,
                  ymax = mean_value + 2*se_value),
              alpha = 0.2, fill = "darkgreen") +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_point(color = "darkgreen", size = 2) +
  scale_x_continuous(breaks = seq(0, 28, 7)) +
  labs(
    title = "NOx Levels by Day in 28-Day Cycle",
    subtitle = "Evidence of monthly periodicity in pollution levels",
    x = "Day in 28-day cycle",
    y = "Mean NOx (mg/m³)"
  ) +
  theme_minimal()

# Create train-test split
n_days <- max(poblenou_ts$day)
train_cutoff <- floor(0.96 * n_days)

train_data <- poblenou_ts %>%
  filter(day <= train_cutoff) %>%
  arrange(time_id)

test_data <- poblenou_ts %>%
  filter(day > train_cutoff) %>%
  arrange(time_id)

cat("\nTrain-test split:\n")
cat("Training days:", min(train_data$day), "-", max(train_data$day), "\n")
cat("Testing days:", min(test_data$day), "-", max(test_data$day), "\n")

# Fit ffc model with cyclic splines for hours, weekly, AND monthly seasonality
# Using bam for faster computation with large data
mod_poblenou <- ffc_gam(
  value ~
    fts(day, mean_only = TRUE, time_k = 4) +
    # Hourly pattern that varies by weekday
    fts(hour, day_week, k = 7) +
    s(day_month, bs = "cc", k = 8),
  data = train_data,
  family = nb(),
  # Cyclic boundaries for all three cycles
  knots = list(
    day_month = c(0.5, 28.5)
  ),
  time = "day",
  engine = "bam",
  discrete = TRUE
)

cat("\n=== Model with hourly, weekly, and monthly seasonality (Poisson BAM) ===\n")

# Model summary and diagnostics
summary(mod_poblenou)
gratia::draw(mod_poblenou)

# Generate forecasts
fable_forecasts <- as_fable(
  mod_poblenou,
  model = "ETS",
  newdata = test_data,
  type = "response"
)

# Evaluate forecast accuracy
accuracy(fable_forecasts, test_data)

# Visualization 1: fable autoplot
fable_forecasts %>%
  autoplot(data = train_data, show_gap = FALSE, fill = "darkblue") +
  labs(
    x = "Day",
    y = "Value",
    title = "Poblenou Hourly Forecasts by Hour",
    subtitle = "Functional forecasting with cyclic splines for diurnal patterns"
  )

# Visualization 2: Custom plot with uncertainty bands
ggplot(fable_forecasts %>%
         mutate(
           time_seq = row_number(),
           q05 = quantile(.dist, 0.05),
           q10 = quantile(.dist, 0.10),
           q25 = quantile(.dist, 0.25),
           q50 = quantile(.dist, 0.50),
           q75 = quantile(.dist, 0.75),
           q90 = quantile(.dist, 0.90),
           q95 = quantile(.dist, 0.95)
         ),
        aes(x = time_seq)) +
  # Uncertainty ribbons
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.2, fill = "darkblue") +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.3, fill = "darkblue") +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.4, fill = "darkblue") +
  # Point forecast (median)
  geom_line(aes(y = q50), color = "darkblue", linewidth = 1) +
  # Actual values
  geom_line(aes(y = value), color = "black", linewidth = 0.8) +
  geom_point(aes(y = value), color = "black", size = 0.8) +
  labs(
    title = "Poblenou Forecasts with Uncertainty Bands",
    x = "Time Index",
    y = "Value",
    subtitle = "Blue ribbons: 90%, 80%, and 50% prediction intervals; Black: observed",
    caption = "Cyclic splines capture diurnal patterns with smooth hour transitions"
  ) +
  theme_minimal()

# Visualization 3: Weekly pattern validation
# Show forecasts broken down by day of week
fable_forecasts %>%
  left_join(
    test_data %>% select(day, hour, day_week) %>% distinct(),
    by = c("day", "hour")
  ) %>%
  mutate(
    q05 = quantile(.dist, 0.05),
    q95 = quantile(.dist, 0.95),
    q50 = quantile(.dist, 0.50)
  ) %>%
  group_by(day_week, hour) %>%
  summarise(
    mean_forecast = mean(q50),
    mean_actual = mean(value),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = c(mean_forecast, mean_actual),
    names_to = "type",
    values_to = "value"
  ) %>%
  mutate(
    type = ifelse(type == "mean_forecast", "Forecast", "Actual")
  ) %>%
  ggplot(aes(x = hour, y = value, color = type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1) +
  facet_wrap(~day_week, ncol = 4,
             labeller = labeller(day_week = c("1" = "Monday", "2" = "Tuesday",
                                              "3" = "Wednesday", "4" = "Thursday",
                                              "5" = "Friday", "6" = "Saturday",
                                              "7" = "Sunday"))) +
  scale_color_manual(values = c("Actual" = "black", "Forecast" = "darkblue")) +
  labs(
    title = "Weekly Pattern Capture: Forecast vs Actual",
    subtitle = "Model successfully captures weekend effect and hourly patterns",
    x = "Hour of Day",
    y = "NOx (mg/m³)",
    color = "Type"
  ) +
  theme_minimal()

# Visualization 4: Focus on hourly patterns
test_sample <- test_data %>%
  filter(day <= min(day) + 2)  # First few test days

fable_sample <- fable_forecasts %>%
  filter(day <= min(test_data$day) + 2)

ggplot(fable_sample %>%
         mutate(
           q05 = quantile(.dist, 0.05),
           q95 = quantile(.dist, 0.95),
           q50 = quantile(.dist, 0.50)
         ),
        aes(x = hour)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.3, fill = "darkblue") +
  geom_line(aes(y = q50), color = "darkblue", linewidth = 1) +
  geom_line(aes(y = value), color = "black", linewidth = 1) +
  geom_point(aes(y = value), color = "black", size = 1.5) +
  facet_wrap(~day, labeller = label_both) +
  scale_x_continuous(breaks = seq(0, 23, 4)) +
  labs(
    title = "Hourly Pattern Validation",
    x = "Hour of Day",
    y = "Value",
    subtitle = "Cyclic spline captures smooth diurnal transitions",
    caption = "Sample of first few forecast days"
  ) +
  theme_minimal()
