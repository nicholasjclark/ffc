# Forecasting from a ffc_gam object
library(ffc)
library(fable)
library(ggplot2); theme_set(theme_classic())

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
 transform <- posreal; gam_fam <- tw()
# transform <- count; gam_fam <- nb()
# transform <- proportional; gam_fam <- betar()
# transform <- binary; gam_fam <- binomial()

# Simulated training and testing data; use Student-T
# for more realism
simdat <- mvgam::sim_mvgam(
  n_series = 1,
  trend_model = mvgam::GP(),
  drift = TRUE,
  prop_trend = 0.6,
  prop_train = 0.7,
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
  dplyr::filter(time <= 75)

test_tsibble <- dat_tsibble %>%
  dplyr::filter(time > 75)

# Fit a ffc model
mod <- ffc_gam(
  y ~
    # Use mean_only = TRUE to ensure that only a constant
    # basis function is used as the by-variable to fit a smooth
    # of 'time' that we can then forecast ahead
    fts(
      time,
      mean_only = TRUE,
      time_k = 60
    ) +
    fts(
      season,
      bs = 'cc',
      k = 8
    ),
  knots = list(season = c(0.5, 12.5)),
  time = "time",
  data = train_tsibble,
  family = gam_fam
)
summary(mod)

# View the time-varying basis coefficients
gratia::draw(mod)

# Compute forecast distribution by fitting the basis coefficient
# time series models

# First try GP factors
fcgp <- forecast(
  object = mod,
  newdata = test_tsibble,
  model = 'GPDF',
  K = 2,
  # use summary = FALSE to return the full distribution
  summary = FALSE
)

# Convert resulting forecasts to a fable for automatic
# plotting / scoring of forecasts
newdata <- test_tsibble
newdata[['value']] <- fcgp

fc_ffcgp <- fabletools:::build_fable(
  newdata,
  response = 'y',
  distribution = 'value'
)

# Now try with AR(4) factors
fcar <- forecast(
  object = mod,
  newdata = test_tsibble,
  model = 'ARDF',
  K = 2,
  lag = 4, # AR order
  summary = FALSE
)
newdata[['value']] <- fcar
fc_ffcar <- fabletools:::build_fable(
  newdata,
  response = 'y',
  distribution = 'value'
)

# Now try with VAR(2) factors
fcvar <- forecast(
  object = mod,
  newdata = test_tsibble,
  model = 'VARDF',
  K = 3, # number of factors
  lag = 2, # AR order
  summary = FALSE
)
newdata[['value']] <- fcvar
fc_ffcvar <- fabletools:::build_fable(
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
fc_ffcgp %>%
  autoplot(train_tsibble) +
  geom_line(
    data = test_tsibble,
    aes(y = y)
  ) +
  ggtitle('FFC GP forecast')

fc_ffcar %>%
  autoplot(train_tsibble) +
  geom_line(
    data = test_tsibble,
    aes(y = y)
  ) +
  ggtitle('FFC AR forecast')

fc_ffcvar %>%
  autoplot(train_tsibble) +
  geom_line(
    data = test_tsibble,
    aes(y = y)
  ) +
  ggtitle('FFC VAR forecast')

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
fc_ffcgp %>%
  fabletools::accuracy(
    test_tsibble,
    measures = distribution_accuracy_measures
  )

fc_ffcar %>%
  fabletools::accuracy(
    test_tsibble,
    measures = distribution_accuracy_measures
  )

fc_ffcvar %>%
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


# A spatiotemporal species distribution example
# Download the pcod data from the sdmTMB package
temp <- tempfile()
download.file('https://github.com/cran/sdmTMB/raw/refs/heads/master/data/pcod.rda',
              temp)
load(temp)
unlink(temp)

unique(pcod$year)
data_train <- pcod %>%
  dplyr::filter(year < 2015)
data_test <- pcod %>%
  dplyr::filter(year >= 2015)

# Model pcod presence / absence over space and time
mod <- ffc_gam(
  present ~
    s(depth_scaled, k = 3) +
    fts(lon, lat, mean_only = TRUE,
        time_k = 7) +
    fts(lon, lat, k = 6,
        time_k = 7),
  data = data_train,
  time = 'year',
  family = binomial(),
  engine = 'gam',
  method = 'REML'
)
summary(mod)
gratia::draw(mod)

# Forecast and calculate probability of occurrence (expectation)
fc <- forecast(
  object = mod,
  newdata = data_test,
  model = 'GPDF',
  K = 3,
  type = 'expected',
  summary = TRUE
)

# Plot true out of sample occurrences
ggplot(data_test,
       aes(x = lat,
           y = lon,
           colour = as.factor(present))) +
  facet_wrap(~year) +
  geom_point() +
  scale_colour_viridis_d() +
  theme_bw()

# Plot predicted probability of occurrence
ggplot(fc %>%
         dplyr::bind_cols(data_test),
       aes(x = lat,
           y = lon,
           colour = .estimate)) +
  facet_wrap(~year) +
  geom_point() +
  scale_colour_viridis_c() +
  theme_bw()


# Now try modelling density over space and time
mod <- ffc_gam(
  density ~
    s(depth_scaled, k = 3) +
    fts(lon, lat,
        mean_only = TRUE,
        time_k = 7) +
    fts(lon, lat,
        k = 6,
        time_k = 7),
  data = data_train,
  time = 'year',
  family = tw(),
  engine = 'gam',
  method = 'REML'
)
summary(mod)
gratia::draw(mod)

# Forecast
fc <- forecast(
  object = mod,
  newdata = data_test,
  model = 'GPDF',
  K = 3,
  type = 'response',
  summary = TRUE
)

ggplot(data_test,
       aes(x = lat,
           y = lon,
           colour = log(density + 1))) +
  facet_wrap(~year) +
  geom_point() +
  scale_colour_viridis_c() +
  theme_bw()

ggplot(fc %>%
         dplyr::bind_cols(data_test),
       aes(x = lat,
           y = lon,
           colour = log(.estimate + 1))) +
  facet_wrap(~year) +
  geom_point() +
  scale_colour_viridis_c() +
  theme_bw()

# A seasonal forecasting example
# Load the AirPassengers dataset and plot the time series
data("AirPassengers")
plot(AirPassengers,
     bty = "l",
     lwd = 2,
     col = "darkred"
)

# This plot suggests that the seasonal pattern changes over time,
# not just in magnitude but also perhaps a bit in shape. Convert
# to a data.frame() object
airdat <- mvgam::series_to_mvgam(
  AirPassengers,
  freq = frequency(AirPassengers)
)
dplyr::glimpse(airdat$data_train)

# Plot features of the series
mvgam::plot_mvgam_series(
  data = airdat$data_train,
  newdata = airdat$data_test
)

# Now fit a model that allows both the level and
# the seasonal shape to change over time
mod <- ffc_gam(
  y ~ fts(year, mean_only = TRUE,
          time_k = 25) +
    fts(season, bs = 'cc', k = 12,
        time_k = 10),
  data = airdat$data_train,
  family = nb(),
  time = 'time'
)
summary(mod)
gratia::draw(mod)

# Forecast
fc <- forecast(
  object = mod,
  newdata = airdat$data_test,
  model = 'GPDF',
  summary = TRUE
)

# Plot forecasts
plot_dat <- airdat$data_train %>%
  dplyr::bind_rows(airdat$data_test %>%
                     dplyr::bind_cols(fc))

ggplot(
  plot_dat,
  aes(x = time,
      y = y)
) +
  geom_ribbon(aes(ymax = .q97.5,
                  ymin = .q2.5),
              alpha = 0.15) +
  geom_ribbon(aes(ymax = .q90,
                  ymin = .q10),
              alpha = 0.2) +
  geom_line()
