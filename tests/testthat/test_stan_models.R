# Tests for Stan-based dynamic factor models (ARDF, GPDF, VARDF)
# These tests are skipped on CRAN due to computational intensity

test_that("ARDF model works with basic functionality", {
  skip_on_cran()
  skip_if_not_installed("fable")
  
  # Simple test data
  test_data <- data.frame(
    y = c(10, 12, 14, 13, 15, 16, 18, 17, 19, 20, 
          22, 21, 23, 24, 26, 25, 27, 28, 30, 29),
    season = rep(1:4, 5),
    time = 1:20
  )
  
  # Fit simple model
  mod <- SW(ffc_gam(
    y ~ fts(time, mean_only = TRUE, time_k = 8, k = 3),
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
    iter = 200,
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
  skip_if_not_installed("fable")
  
  # Simple test data
  test_data <- data.frame(
    y = c(5, 7, 9, 8, 10, 11, 13, 12, 14, 15,
          17, 16, 18, 19, 21, 20, 22, 23, 25, 24),
    season = rep(1:4, 5),
    time = 1:20
  )
  
  # Fit simple model
  mod <- SW(ffc_gam(
    y ~ fts(time, mean_only = TRUE, time_k = 6, k = 2),
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
  skip_if_not_installed("fable")
  
  # Simple test data with more observations for VAR stability
  test_data <- data.frame(
    y = c(8, 10, 12, 11, 13, 14, 16, 15, 17, 18,
          20, 19, 21, 22, 24, 23, 25, 26, 28, 27,
          29, 30, 32, 31, 33, 34, 36, 35, 37, 38),
    season = rep(1:4, 7.5)[1:30],
    time = 1:30
  )
  
  # Fit model with multiple basis functions for VAR
  mod <- SW(ffc_gam(
    y ~ fts(time, mean_only = TRUE, time_k = 10, k = 3),
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

test_that("Stan models work with different families", {
  skip_on_cran()
  skip_if_not_installed("fable")
  
  # Test with Poisson family (count data)
  test_data <- data.frame(
    y = rpois(15, lambda = exp(1 + 0.1 * (1:15))),
    time = 1:15
  )
  
  mod_poisson <- SW(ffc_gam(
    y ~ fts(time, mean_only = TRUE, time_k = 6, k = 2),
    time = "time",
    data = test_data,
    family = poisson()
  ))
  
  newdata <- data.frame(time = c(16, 17))
  
  # Test GPDF with Poisson
  fc_poisson <- SW(forecast(
    mod_poisson,
    newdata = newdata,
    model = "GPDF",
    K = 1,
    chains = 1, 
    iter = 200,
    summary = TRUE
  ))
  
  expect_true(inherits(fc_poisson, "tbl_df"))
  expect_true(all(fc_poisson$.estimate >= 0))  # Poisson should give non-negative predictions
})

test_that("Stan models validate parameters correctly", {
  skip_on_cran()
  
  # Simple model for parameter testing
  test_data <- data.frame(
    y = 1:10,
    time = 1:10  
  )
  
  mod <- SW(ffc_gam(
    y ~ fts(time, mean_only = TRUE, time_k = 5, k = 2),
    time = "time",
    data = test_data,
    family = gaussian()
  ))
  
  newdata <- data.frame(time = 11)
  
  # Test invalid K parameter
  expect_error(
    forecast(mod, newdata = newdata, model = "ARDF", K = 0),
    "K must be positive"
  )
  
  # Test invalid lag parameter  
  expect_error(
    forecast(mod, newdata = newdata, model = "ARDF", K = 1, lag = 0),
    "lag must be positive"
  )
})

test_that("Stan models return distributional objects when summary=FALSE", {
  skip_on_cran()
  skip_if_not_installed("fable")
  skip_if_not_installed("distributional")
  
  test_data <- data.frame(
    y = c(5, 7, 9, 11, 13, 15),
    time = 1:6
  )
  
  mod <- SW(ffc_gam(
    y ~ fts(time, mean_only = TRUE, time_k = 4, k = 2),
    time = "time",
    data = test_data,
    family = gaussian()
  ))
  
  newdata <- data.frame(time = 7)
  
  # Test distributional return
  fc_dist <- SW(forecast(
    mod,
    newdata = newdata,
    model = "GPDF", 
    K = 1,
    chains = 1,
    iter = 200,
    summary = FALSE
  ))
  
  expect_true(inherits(fc_dist, "distribution"))
  expect_equal(length(fc_dist), nrow(newdata))
})

test_that("as_fable works with Stan model forecasts", {
  skip_on_cran()
  skip_if_not_installed("fable")
  skip_if_not_installed("distributional")
  
  test_data <- data.frame(
    y = c(2, 4, 6, 8, 10, 12),
    time = 1:6
  )
  
  mod <- SW(ffc_gam(
    y ~ fts(time, mean_only = TRUE, time_k = 4, k = 2), 
    time = "time",
    data = test_data,
    family = gaussian()
  ))
  
  newdata <- data.frame(
    y = c(14, 16),
    time = c(7, 8)
  )
  
  # Generate forecasts with ARDF
  forecasts <- SW(forecast(
    mod,
    newdata = newdata,
    model = "ARDF",
    K = 1,
    lag = 1,
    chains = 1,
    iter = 200,
    summary = FALSE
  ))
  
  # Convert to fable
  fc_fable <- SW(as_fable(mod, newdata = newdata, forecasts = forecasts))
  
  expect_true(inherits(fc_fable, "fbl_ts"))
  expect_true(inherits(fc_fable$.dist, "distribution")) 
  expect_equal(nrow(fc_fable), nrow(newdata))
})