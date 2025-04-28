# Forecasting from a ffc_gam object
library(ffc)
library(ggplot2); theme_set(theme_classic())

# Simulated training and testing data
set.seed(0)
simdat <- mvgam::sim_mvgam(
  n_series = 1,
  trend_model = mvgam::GP(),
  prop_trend = 0.6,
  prop_train = 0.85,
  mu = 2,
  family = nb()
)
mvgam::plot_mvgam_series(
  data = simdat$data_train,
  newdata = simdat$data_test
)

# Fit a model
mod <- ffc_gam(
  y ~ s(season, bs = "cc", k = 12) +

    # Use mean_only = TRUE to ensure that only a constant
    # basis function is used as the by-variable; this enables
    # us to fit a smooth of 'time' that we can then forecast
    # ahead just like any other fts basis :)
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
gratia::draw(mod)
autoplot(
  fts_coefs(
    mod,
    summary = FALSE,
    times = 10
  )
)

# Compute forecast distribution
fc <- forecast(
  mod,
  newdata = simdat$data_test,
  n_sims = 100
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
    alpha = 0.4,
    fill = 'darkred'
  ) +
  geom_line(
    aes(y = .estimate),
    col = 'darkred',
    linewidth = 1
  ) +
  geom_point(
    pch = 21,
    fill = 'black',
    col = 'white'
  )


# Compare to a model that uses spline extrapolation instead
mod2 <- ffc_gam(
  y ~ s(season, bs = "cc", k = 12) +
    s(time, k = 30),

  # Supply some knots to ensure they are picked up correctly
  knots = list(season = c(0.5, 12.5)),
  time = "time",
  data = simdat$data_train,
  family = nb()
)
summary(mod2)
gratia::draw(mod2)

# Compute forecast distribution
fc2 <- forecast(
  mod2,
  newdata = simdat$data_test,
  n_sims = 100
)

# Plot forecasts against truth
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
    alpha = 0.4,
    fill = 'darkred'
  ) +
  geom_line(
    aes(y = .estimate),
    col = 'darkred',
    linewidth = 1
  ) +
  geom_point(
    pch = 21,
    fill = 'black',
    col = 'white'
  )
