# Setup models for tests
library("testthat")
library("ffc")
data("growth_data")

# Fit some models that can be used for testing several downstream functions
example_mod <- ffc_gam(
  formula = height_cm ~
    s(id, bs = "re") +
    fts(
      age_yr,
      bs = "tp", k = 4,
      time_bs = "cr", time_k = 6
    ),
  data = growth_data,
  time = "age_yr",
  family = mgcv::Tweedie(p = 1.5),
  engine = "bam",
  discrete = TRUE
)

# Another example, this time for a time series that only aims to forecast
# the temporal trend smooth
dat <- structure(
  list(
    y = c(
      38L, 35L, 26L, 31L, 67L, 58L, 42L, 28L, 20L, 29L, 38L, 22L, 18L,
      12L, 8L, 15L, 13L, 12L, 9L, 7L, 8L,
      13L, 16L, 15L, 13L, 18L, 17L, 20L, 40L, 50L, 26L, 12L, 15L, 46L,
      52L, 26L, 13L, 21L, 14L, 10L, 24L, 20L, 15L, 6L, 4L, 7L, 13L,
      13L, 4L, 3L, 8L, 11L, 21L, 15L, 11L, 9L, 3L, 24L, 20L, 4L, 7L,
      4L, 11L, 15L, 27L, 25L, 14L, 10L, 14L, 17L, 16L, 13L, 11L, 7L,
      16L
    ),
    season = c(
      1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L,
      12L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 1L, 2L,
      3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 1L, 2L, 3L, 4L, 5L,
      6L, 7L, 8L, 9L, 10L, 11L, 12L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L,
      9L, 10L, 11L, 12L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L,
      12L, 1L, 2L, 3L
    ),
    year = c(
      1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
      1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L,
      3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L,
      4L, 4L, 4L, 4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L,
      5L, 5L, 5L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 7L,
      7L, 7L
    ),
    time = 1:75
  ),
  row.names = c(NA, -75L),
  class = "data.frame"
)

example_mod2 <- ffc_gam(
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
  data = dat,
  family = poisson()
)

# Some utility testing functions
expect_match2 <- function(object, regexp) {
  any(grepl(regexp, object, fixed = TRUE))
}

expect_character <- function(object, ...) {
  testthat::expect_true(is(object, "character"), ...)
}

expect_list <- function(object, ...) {
  testthat::expect_true(is(object, "list"), ...)
}

expect_ggplot <- function(object, ...) {
  testthat::expect_true(is(object, "ggplot"), ...)
}

expect_range <- function(object, lower = -Inf, upper = Inf, ...) {
  testthat::expect_true(all(object >= lower & object <= upper), ...)
}

SM <- suppressMessages
SW <- suppressWarnings
