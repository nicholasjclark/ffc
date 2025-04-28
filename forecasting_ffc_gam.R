# Forecasting from a ffc_gam object
library(ffc)
library(ggplot2); theme_set(theme_classic())

# Simulated training and testing data
set.seed(0)
simdat <- mvgam::sim_mvgam(
  n_series = 1,
  trend_model = mvgam::GP(),
  prop_trend = 0.65,
  prop_train = 0.80,
  mu = 1.5,
  family = nb()
)
mvgam::plot_mvgam_series(
  data = simdat$data_train,
  newdata = simdat$data_test
)

# Fit a model
mod <- ffc_gam(
  y ~
    # Capture possible time-varying seasonality
    fts(
      season,
      bs = "cc",
      k = 8
    ) +

    # Use mean_only = TRUE to ensure that only a constant
    # basis function is used as the by-variable to fit a smooth
    # of 'time' that we can then forecast ahead
    fts(
      time,
      mean_only = TRUE,
      time_k = 30
    ),

  # Supply some knots to ensure they are picked up correctly
  knots = list(season = c(0.5, 12.5)),
  time = "time",
  data = simdat$data_train,
  family = nb()
)
summary(mod)

# View draws of the time-varying basis coefficient
autoplot(
  fts_coefs(
    mod,
    summary = FALSE,
    times = 10
  )
)

# Compute forecast distribution by fitting the basis coefficient
# time series models in parallel (which is automatically supported
# within the fable package)
library(future)
plan(multisession)

fc <- forecast(
  mod,
  newdata = simdat$data_test,
  stationary = TRUE
)

# Plot forecasts against truth
plot_dat <- dplyr::bind_rows(
  simdat$data_train,
  (simdat$data_test %>%
     dplyr::bind_cols(fc))
)

ggplot(
  data = plot_dat,
  aes(x = time,
      y = y)
) +
  geom_ribbon(
    aes(ymax = .q97.5,
        ymin = .q2.5),
    alpha = 0.3,
    fill = 'darkred'
  ) +
  geom_ribbon(
    aes(ymax = .q90,
        ymin = .q10),
    alpha = 0.4,
    fill = 'darkred'
  ) +
  geom_line(
    aes(y = .estimate),
    col = 'darkred',
    linewidth = 1
  ) +
  geom_line()

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

sim_tsibble %>%
  model(
    arima = ARIMA(fabletools::box_cox(y, feasts::guerrero(y)))
  ) %>%
  forecast(
    h = max(simdat$data_test$time) -
      min(simdat$data_test$time) + 1
  ) %>%
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
  n_sims = 100
)

plot_dat <- dplyr::bind_rows(
  simdat$data_train,
  (simdat$data_test %>%
     dplyr::bind_cols(fc2))
)

ggplot(
  data = plot_dat,
  aes(x = time,
      y = y)
) +
  geom_ribbon(
    aes(ymax = .q97.5,
        ymin = .q2.5),
    alpha = 0.3,
    fill = 'darkred'
  ) +
  geom_ribbon(
    aes(ymax = .q90,
        ymin = .q10),
    alpha = 0.4,
    fill = 'darkred'
  ) +
  geom_line(
    aes(y = .estimate),
    col = 'darkred',
    linewidth = 1
  ) +
  geom_line()
