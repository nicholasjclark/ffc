# Setup models for tests
library("testthat")
library("ffc")
data("growth_data")

# Fit a model that can be used for testing several downstream functions
example_mod <- ffc_gam(
  formula = height_cm ~
    s(id, bs = 're') +
    fts(age_yr, bs = 'tp', k = 4,
        time_bs = 'cr', time_k = 6),
  data = growth_data,
  time = "age_yr",
  family = mgcv::Tweedie(p = 1.5),
  engine = 'bam',
  discrete = TRUE
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
