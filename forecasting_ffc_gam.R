# Forecasting from a ffc_gam object
library(ffc)
simdat <- mvgam::sim_mvgam(
  n_series = 1,
  trend_model = mvgam::RW(),
  prop_trend = 0.75,
  mu = 2
)
mod <- ffc_gam(
  y ~ s(season, bs = "cc", k = 12) +

    # Use mean_only = TRUE to ensure that only a constant
    # basis function is used as the by-variable; this enables
    # us to fit a smooth of 'time' that we can then forecast
    # ahead just like any other fts basis :)
    fts(
      time,
      mean_only = TRUE,
      time_k = 20, time_bs = "bs", time_m = 1
    ),

  # Supply some knots to ensure they are picked up correctly
  knots = list(season = c(0.5, 12.5)),
  time = "time",
  data = simdat$data_train,
  family = poisson()
)

functional_coefs <- fts_coefs(
  mod,
  summary = FALSE,
  times = 10
)

functional_fc <- forecast(
  object = functional_coefs,
  h = 5,
  times = 5,
  model = ETS()
)
functional_fc

fc <- forecast(
  mod,
  newdata = simdat$data_test
)
