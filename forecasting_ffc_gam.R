# Forecasting from a ffc_gam object
library(ffc)
library(fable)
library(ggplot2); theme_set(theme_classic())
library(future)
plan(multisession, workers = 4)

# Define transformations that avoid distributional
# assumptions so that ffc is not unfairly advantaged
# over ARIMA / Box-Cox
count = function(x) {
  xhat <- x * runif(1, 1, 20)
  floor(xhat + abs(min(xhat)))
}

proportional = function(x) {
  plogis(x)
}

posreal = function(x) {
  xhat <- x * runif(1, 1, 200)
  xhat + abs(min(xhat)) + runif(1, 0.01, 5)
}

binary = function(x) {
  rbinom(
    n = length(x),
    size = 1,
    prob = plogis(x)
  )
}

# Simulation setup
# transform <- posreal; gam_fam <- tw()
transform <- count; gam_fam <- nb()
# transform <- proportional; gam_fam <- betar()
# transform <- binary; gam_fam <- binomial()

# Simulated training and testing data; use Student-T
# for more realism
simdat <- mvgam::sim_mvgam(
  n_series = 1,
  trend_model = mvgam::RW(),
  drift = TRUE,
  prop_trend = 0.6,
  prop_train = 0.9,
  mu = -1,
  family = mvgam::student()
)
mvgam::plot_mvgam_series(
  data = simdat$data_train,
  newdata = simdat$data_test
)

# Convert to non-Gaussian
dat_tsibble <- rbind(
  simdat$data_train,
  simdat$data_test
) %>%
  dplyr::mutate(
    yearmon = tsibble::make_yearmonth(
      year = year,
      month = season
    )
  ) %>%
  tsibble::as_tsibble(
    index = yearmon
  ) %>%
  dplyr::mutate(y = do.call(
    transform,
    list(y)
  )
  )

# Create tsibbles of training and testing data
train_tsibble <- dat_tsibble %>%
  dplyr::filter(time <= 85)

test_tsibble <- dat_tsibble %>%
  dplyr::filter(time > 85)

# Fit a ffc model
mod <- ffc_gam(
  y ~
    # Use mean_only = TRUE to ensure that only a constant
    # basis function is used as the by-variable to fit a smooth
    # of 'time' that we can then forecast ahead
    fts(
      time,
      mean_only = TRUE,
      time_bs = 'tp',
      time_k = 50
    ) +
    fts(
      season,
      bs = 'cc',
      k = 6,
      time_bs = 'tp'
    ),
  knots = list(season = c(0.5, 12.5)),
  time = "time",
  data = train_tsibble,
  family = gam_fam,
  select = TRUE
)
summary(mod)

# View the time-varying basis coefficients
gratia::draw(mod)

# Compute forecast distribution by fitting the basis coefficient
# time series models in parallel (which is automatically supported
# within the fable package)
fc <- forecast(
  object = mod,
  newdata = test_tsibble,
  model = 'ETS',
  quantile_fc = FALSE,
  stationary = FALSE,
  # use summary = FALSE to return the full distribution
  summary = FALSE
)

# Convert resulting forecasts to a fable for automatic
# plotting / scoring of forecasts
newdata <- test_tsibble
newdata[['value']] <- fc

fc_ffc <- fabletools:::build_fable(
  newdata,
  response = 'y',
  distribution = 'value'
)

# Benchmarks
fc_arima <- train_tsibble %>%
  model(
    arima = ARIMA(fabletools::box_cox(y, feasts::guerrero(y)))
  ) %>%
  forecast(
    h = max(test_tsibble$time) -
      min(test_tsibble$time) + 1
  )

fc_ets <- train_tsibble %>%
  model(
    ets = ETS(fabletools::box_cox(y, feasts::guerrero(y)))
  ) %>%
  forecast(
    h = max(test_tsibble$time) -
      min(test_tsibble$time) + 1
  )

# Plot the resulting forecasts
fc_ffc %>%
  autoplot(train_tsibble) +
  geom_line(
    data = test_tsibble,
    aes(y = y)
  ) +
  ggtitle('FFC forecast')

fc_arima %>%
  autoplot(train_tsibble) +
  geom_line(
    data = test_tsibble,
    aes(y = y)
  ) +
  ggtitle('ARIMA with Box-Cox forecast')

fc_ets %>%
  autoplot(train_tsibble) +
  geom_line(
    data = test_tsibble,
    aes(y = y)
  ) +
  ggtitle('ETS with Box-Cox forecast')

# Score the forecasts with CRPS (lower is better)
fc_ffc %>%
  fabletools::accuracy(
    test_tsibble,
    measures = distribution_accuracy_measures
  )

fc_arima %>%
  fabletools::accuracy(
    test_tsibble,
    measures = distribution_accuracy_measures
  )

fc_ets %>%
  fabletools::accuracy(
    test_tsibble,
    measures = distribution_accuracy_measures
  )
