library(mvgam)
library(tidyverse)
library(fable)
library(tsibble)

# Get a multivariate time series
series_dat <- tourism %>%
  filter(
    Region == "Melbourne"
  ) %>%
  mutate(
    quarter = lubridate::quarter(Quarter)
  ) %>%
  # convert to tibble for easier selecting of variables
  as_tibble() %>%
  group_by(Purpose) %>%
  mutate(
    time = 1:n()
  ) %>%
  ungroup() %>%
  # add 'series' indicator
  mutate(series = as.factor(Purpose),
         y = Trips) %>%
  arrange(series, time) %>%
  # series_dat should already have the training and
  # testing times, with NAs for the out of sample obs
  mutate(y = case_when(
    time <= 70 ~ y,
    TRUE ~ NA
  )) %>%
  select(y, series, time)

# Arguments
n_lv <- 2
family <- 1
prior_alpha <- c(mean(series_dat$y,
                      na.rm = TRUE),
                 50)
prior_sigma <- c(4, 1, 4)
prior_ar <- c(1, 1, 1)
p <- 2
series_dat <- series_dat

# Pull some key arguments
n_series <- nlevels(series_dat$series)
n_timepoints <- length(unique(series_dat$time))

# Need a matrix of values for the series, excluding forecast values
Y <- series_dat %>%
  pivot_wider(names_from = 'series',
              values_from = 'y') %>%
  select(-time) %>%
  as.matrix()

stopifnot(identical(NROW(Y), n_timepoints))
stopifnot(identical(NCOL(Y), n_series))

# Create matrix representing whether an observation was missing or not
Y_observed <- matrix(
  NA,
  ncol = NCOL(Y),
  nrow = NROW(Y)
)
for (i in 1:dim(Y)[1]) {
  for (s in 1:dim(Y)[2]) {
    if (is.na(Y[i, s])) {
      Y_observed[i, s] = 0
    } else {
      Y_observed[i, s] = 1
    }
  }
}

# Use -1 for any missing observations so Stan doesn't throw errors due to NAs
Y[is.na(Y)] <- -1

# model data list
model_data <- list(
  n = n_timepoints,
  n_lv = n_lv,
  n_series = n_series,

  # Will need to compute this as well
  sample_var = c(0.3, 0.3, 0.9, 2.9),
  M = n_lv * (n_series - n_lv) + n_lv * (n_lv - 1) / 2 + n_lv,
  n_nonmissing = length(which(Y_observed == 1)),
  flat_ys = as.vector(Y)[which(
    as.vector(Y_observed) == 1
  )],
  obs_ind = which(as.vector(Y_observed) == 1),
  family = family,
  prior_alpha = prior_alpha,
  prior_sigma = prior_sigma,
  prior_ar = prior_ar,
  P = p,
  beta = array(1)
)

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
stan_mod <- rstan::stan_model('inst/stan/ar_df.stan')

fit1 <- rstan::sampling(
  stan_mod,
  data = model_data,
  chains = 4,
  cores = 4,
  refresh = 100,
  iter = 800,
  warmup = 400
)

out <- summary(
  fit1,
  pars = c('alpha', 'sigma_obs', 'ar', 'nu', 'Lambda')
)
out$summary

# Extract predictions for each series
preds <- as.matrix(fit1, pars = 'ypred')
ends <- seq(
  0,
  NCOL(preds),
  length.out = n_series + 1
)
starts <- ends + 1
starts <- c(1, starts[-c(1, (n_series + 1))])
ends <- ends[-1]

series_preds <- lapply(
  1:n_series,
  function(x){
    preds[ , starts[x]:ends[x]]
  }
)

plot_draws = function(series_preds, series = 1){
  s_preds <- series_preds[[series]]
  creds <- apply(s_preds, 2, function(x)
    quantile(x, probs = c(0.025, 0.5, 0.975)))
  s_name <- levels(series_dat$series)[series]
  ylim <- range(s_preds)
  plot(
    s_preds[1, ],
    type = 'l',
    ylim = ylim,
    col = 'white'
  )

  lines(
    creds[2,],
    col = 'grey20',
    lwd = 1.5
  )
  lines(
    creds[1,],
    col = 'grey50',
    lwd = 1.5
  )
  lines(
    creds[3,],
    col = 'grey50',
    lwd = 1.5
  )
  points(
    series_dat %>%
      dplyr::filter(series == s_name) %>%
      dplyr::pull(y),
    pch = 16
  )
}

plot_draws(
  series_preds,
  series = 1
)

plot_draws(
  series_preds,
  series = 2
)

plot_draws(
  series_preds,
  series = 3
)

plot_draws(
  series_preds,
  series = 4
)
