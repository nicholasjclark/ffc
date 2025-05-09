# Forecasting from a ffc_gam object
library(ffc)
library(ggplot2); theme_set(theme_classic())

# Simulated training and testing data
simdat <- mvgam::sim_mvgam(
  n_series = 1,
  trend_model = mvgam::GP(),
  drift = TRUE,
  prop_trend = 0.60,
  prop_train = 0.85,
  mu = 1.5,
  family = poisson()
)
mvgam::plot_mvgam_series(
  data = simdat$data_train,
  newdata = simdat$data_test
)

# Fit a model
mod <- ffc_gam(
  y ~
    # Use mean_only = TRUE to ensure that only a constant
    # basis function is used as the by-variable to fit a smooth
    # of 'time' that we can then forecast ahead
    fts(
      time,
      mean_only = TRUE,
      time_bs = 'ts',
      time_k = 30
    ) +
    # Capture possible time-varying seasonality
    fts(
      season,
      bs = "cc",
      k = 8,
      time_bs = 'ts',
      time_k = 30
    ),

  # Supply some knots to ensure they are picked up correctly
  knots = list(season = c(0.5, 12.5)),
  time = "time",
  data = simdat$data_train,
  family = nb(),
  engine = 'bam',
  discrete = TRUE
)
summary(mod)

# View draws of the time-varying basis coefficient
func_coefs <- fts_coefs(
  mod,
  summary = FALSE,
  times = 20
)
autoplot(
  func_coefs
)

# Compute forecast distribution by fitting the basis coefficient
# time series models in parallel (which is automatically supported
# within the fable package)
library(future)
plan(multisession)

fc <- forecast(
  object = mod,
  newdata = simdat$data_test,
  model = 'ARIMA',
  stationary = FALSE,
  n_draws = 20,
  n_sims = 100,
  # use summary = FALSE to return the full distribution
  summary = FALSE
)

# Convert resulting forecasts to a fable for automatic
# plotting / scoring of forecasts
newdata <- simdat$data_test %>%
  dplyr::mutate(
    yearmon = tsibble::make_yearmonth(
      year = year,
      month = season
    )
  ) %>%
  tsibble::as_tsibble(
    index = yearmon
  )
newdata[['value']] <- fc

fc_ffc <- fabletools:::build_fable(
  newdata,
  response = 'y',
  distribution = 'value'
)

# How would an automatic ARIMA model with Box-Cox transformation compare?
library(fable)
sim_tsibble <- simdat$data_train %>%
  dplyr::mutate(
    yearmon = tsibble::make_yearmonth(
      year = year,
      month = season
    )
  ) %>%
  tsibble::as_tsibble(
    index = yearmon
  )

fc_arima <- sim_tsibble %>%
  model(
    arima = ARIMA(fabletools::box_cox(y, feasts::guerrero(y)))
  ) %>%
  forecast(
    h = max(simdat$data_test$time) -
      min(simdat$data_test$time) + 1
  )

# How would a model that uses spline extrapolation compare?
mod2 <- ffc_gam(
  y ~ te(
    season,
    time,
    bs = c("cc", "tp"),
    k = c(4, 10)
  ),

  # Supply some knots to ensure they are picked up correctly
  knots = list(season = c(0.5, 12.5)),
  time = "time",
  data = simdat$data_train,
  family = nb()
)

fc2 <- forecast(
  mod2,
  newdata = simdat$data_test,
  n_sims = 100,
  summary = FALSE
)
newdata[['value']] <- fc2

fc_gam <- fabletools:::build_fable(
  newdata,
  response = 'y',
  distribution = 'value'
)

# Plot the resulting forecasts from all three models
fc_ffc %>%
  autoplot(sim_tsibble) +
  geom_line(
    data = newdata,
    aes(y = y)
  ) +
  ggtitle('FFC forecast')

fc_arima %>%
  autoplot(sim_tsibble) +
  geom_line(
    data = simdat$data_test %>%
      dplyr::mutate(
        yearmon = tsibble::make_yearmonth(
          year = year,
          month = season
        )
      ),
    aes(y = y)
  ) +
  ggtitle('ARIMA with Box-Cox forecast')

fc_gam %>%
  autoplot(sim_tsibble) +
  geom_line(
    data = newdata,
    aes(y = y)
  ) +
  ggtitle('GAM forecast')

# Score the forecasts from all three models with CRPS (lower is better)
fc_ffc %>%
  fabletools::accuracy(
    newdata,
    measures = distribution_accuracy_measures
  )

fc_arima %>%
  fabletools::accuracy(
    newdata,
    measures = distribution_accuracy_measures
  )

fc_gam %>%
  fabletools::accuracy(
    newdata,
    measures = distribution_accuracy_measures
  )
