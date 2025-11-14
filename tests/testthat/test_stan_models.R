# Tests for Stan-based dynamic factor models (ARDF, GPDF, VARDF)
# These tests are skipped on CRAN due to computational intensity

test_that("ARDF model works with basic functionality", {
  skip_on_cran()

  # Simple test data
  test_data <- data.frame(
    y = c(10, 12, 14, 13, 15, 16, 18, 17, 19, 20,
          22, 21, 23, 24, 26, 25, 27, 28, 30, 29),
    season = rep(1:4, 5),
    time = 1:20
  )

  # Fit simple model
  mod <- SW(ffc_gam(
    y ~ fts(time, mean_only = TRUE, time_k = 8, k = 3) +
      fts(season, k = 3, time_k = 8),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  # Test ARDF forecasting
  newdata <- data.frame(
    season = c(1, 2),
    time = c(21, 22)
  )

  # Basic ARDF forecast
  fc_ardf <- SW(forecast(
    mod,
    newdata = newdata,
    model = "ARDF",
    K = 1,
    lag = 1,
    chains = 1,
    iter = 220,
    summary = TRUE
  ))

  expect_true(inherits(fc_ardf, "tbl_df"))
  expect_true(all(c(".estimate", ".error") %in% names(fc_ardf)))
  expect_equal(nrow(fc_ardf), nrow(newdata))
  expect_true(all(is.finite(fc_ardf$.estimate)))
  expect_true(all(fc_ardf$.error >= 0))
})

test_that("GPDF model works with basic functionality", {
  skip_on_cran()

  # Simple test data
  test_data <- data.frame(
    y = c(5, 7, 9, 8, 10, 11, 13, 12, 14, 15,
          17, 16, 18, 19, 21, 20, 22, 23, 25, 24),
    season = rep(1:4, 5),
    time = 1:20
  )

  # Fit simple model
  mod <- SW(ffc_gam(
    y ~ fts(time, mean_only = TRUE, time_k = 6, k = 2) +
      fts(season, k = 3, time_k = 8),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  # Test GPDF forecasting
  newdata <- data.frame(
    season = c(1, 2),
    time = c(21, 22)
  )

  # Basic GPDF forecast
  fc_gpdf <- SW(forecast(
    mod,
    newdata = newdata,
    model = "GPDF",
    K = 1,
    chains = 1,
    iter = 200,
    summary = TRUE
  ))

  expect_true(inherits(fc_gpdf, "tbl_df"))
  expect_true(all(c(".estimate", ".error") %in% names(fc_gpdf)))
  expect_equal(nrow(fc_gpdf), nrow(newdata))
  expect_true(all(is.finite(fc_gpdf$.estimate)))
  expect_true(all(fc_gpdf$.error >= 0))
})

test_that("VARDF model works with basic functionality", {
  skip_on_cran()

  # Simple test data with more observations for VAR stability
  test_data <- data.frame(
    y = c(8, 10, 12, 11, 13, 14, 16, 15, 17, 18,
          20, 19, 21, 22, 24, 23, 25, 26, 28, 27,
          29, 30, 32, 31, 33, 34, 36, 35, 37, 38),
    season = rep(1:4, 9)[1:30],
    time = 1:30
  )

  # Fit model with multiple basis functions for VAR
  mod <- SW(ffc_gam(
    y ~ fts(time, mean_only = TRUE, time_k = 10, k = 3) +
      fts(season, k = 3, time_k = 8),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  # Test VARDF forecasting
  newdata <- data.frame(
    season = c(1, 2),
    time = c(31, 32)
  )

  # Basic VARDF forecast
  fc_vardf <- SW(forecast(
    mod,
    newdata = newdata,
    model = "VARDF",
    K = 2,
    lag = 1,
    chains = 1,
    iter = 200,
    summary = TRUE
  ))

  expect_true(inherits(fc_vardf, "tbl_df"))
  expect_true(all(c(".estimate", ".error") %in% names(fc_vardf)))
  expect_equal(nrow(fc_vardf), nrow(newdata))
  expect_true(all(is.finite(fc_vardf$.estimate)))
  expect_true(all(fc_vardf$.error >= 0))
})

test_that("Stan models validate parameters correctly", {
  skip_on_cran()

  # Simple test data
  test_data <- data.frame(
    y = c(5, 7, 9, 8, 10, 11, 13, 12, 14, 15,
          17, 16, 18, 19, 21, 20, 22, 23, 25, 24),
    season = rep(1:4, 5),
    time = 1:20
  )

  mod <- SW(ffc_gam(
    y ~ fts(time, mean_only = TRUE, time_k = 5, k = 3) +
      fts(season, k = 3, time_k = 8),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  newdata <- data.frame(time = 21:22, season = 1:2)

  # Test invalid K parameter
  expect_error(
    forecast(mod, newdata = newdata, model = "ARDF", K = 0),
    "Assertion on 'K' failed: Must be >= 1"
  )

  # Test invalid lag parameter
  expect_error(
    forecast(mod, newdata = newdata, model = "ARDF", K = 1, lag = 0),
    "Assertion on 'lag' failed: Must be >= 1"
  )
})
