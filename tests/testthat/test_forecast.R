# Test forecasting functionality
test_that("forecast.fts_ts() works with ARIMA model", {
  # Extract functional coefficients for forecasting
  functional_coefs <- fts_coefs(example_mod, summary = FALSE, times = 3)

  # Test basic ARIMA forecasting
  fc_arima <- SW(forecast(functional_coefs,
                          model = "ARIMA",
                          h = 2,
                          times = 3))

  # Check basic structure
  expect_true(inherits(fc_arima, "tbl_df"))
  expect_true(all(c(".basis", ".sim", ".rep") %in% names(fc_arima)))
  expect_true(length(unique(fc_arima$.rep)) == 3)
  expect_true(all(fc_arima$.rep %in% 1:3))

  # Check forecast horizon
  expect_true(length(unique(fc_arima$.time)) == 2)

  # Check that forecasts are numeric
  expect_true(is.numeric(fc_arima$.sim))
})

test_that("forecast.fts_ts() works with different fable models", {
  functional_coefs <- fts_coefs(example_mod, summary = FALSE, times = 2)

  # Test ETS model
  fc_ets <- SW(forecast(functional_coefs,
                        model = "ETS",
                        h = 2,
                        times = 2))
  expect_true(inherits(fc_ets, "tbl_df"))
  expect_true(all(c(".basis", ".sim", ".rep") %in% names(fc_ets)))

  # Test NAIVE model
  fc_naive <- SW(forecast(functional_coefs,
                          model = "NAIVE",
                          h = 2,
                          times = 2))
  expect_true(inherits(fc_naive, "tbl_df"))
  expect_true(all(c(".basis", ".sim", ".rep") %in% names(fc_naive)))

  # Test random walk
  fc_rw <- SW(forecast(functional_coefs,
                       model = "RW",
                       h = 2,
                       times = 2))
  expect_true(inherits(fc_rw, "tbl_df"))
  expect_true(all(c(".basis", ".sim", ".rep") %in% names(fc_rw)))
})

test_that("forecast.fts_ts() validates parameters correctly", {
  functional_coefs <- fts_coefs(example_mod, summary = FALSE, times = 2)

  # Test invalid h parameter
  expect_error(forecast(functional_coefs, h = 0))
  expect_error(forecast(functional_coefs, h = -1))
  expect_error(forecast(functional_coefs, h = 1.5))

  # Test invalid times parameter
  expect_error(forecast(functional_coefs, times = 0))
  expect_error(forecast(functional_coefs, times = -1))
  expect_error(forecast(functional_coefs, times = 1.5))

  # Test stationary parameter with non-ARIMA model
  expect_error(forecast(functional_coefs, model = "ETS", stationary = TRUE))
})

test_that("forecast.fts_ts() handles stationary ARIMA correctly", {
  functional_coefs <- fts_coefs(example_mod, summary = FALSE, times = 2)

  # Test stationary ARIMA
  fc_stationary <- SW(forecast(functional_coefs,
                               model = "ARIMA",
                               stationary = TRUE,
                               h = 2,
                               times = 2))

  expect_true(inherits(fc_stationary, "tbl_df"))
  expect_true(all(c(".basis", ".sim", ".rep") %in% names(fc_stationary)))
})

test_that("forecast.ffc_gam() works with basic newdata", {
  # Create test newdata for forecasting
  newdata <- data.frame(
    season = c(1, 2, 3),
    time = c(76, 77, 78)
  )

  # Test basic forecasting with ARIMA model
  fc <- SW(forecast(example_mod,
                    newdata = newdata,
                    model = "ARIMA",
                    summary = TRUE))

  # Check structure when summary = TRUE
  expect_true(inherits(fc, "tbl_df"))
  expect_true(all(c(".estimate", ".error") %in% names(fc)))
  expect_equal(nrow(fc), nrow(newdata))

  # Check that estimates are numeric
  expect_true(is.numeric(fc$.estimate))
  expect_true(is.numeric(fc$.error))

  # Check that errors are non-negative
  expect_range(fc$.error, lower = 0)
})

test_that("forecast.ffc_gam() works with different prediction types", {
  newdata <- data.frame(
    season = c(1, 2),
    time = c(76, 77)
  )

  # Test different types
  fc_link <- SW(forecast(example_mod,
                         newdata = newdata,
                         type = "link"))
  fc_expected <- SW(forecast(example_mod,
                             newdata = newdata,
                             type = "expected"))
  fc_response <- SW(forecast(example_mod,
                             newdata = newdata,
                             type = "response"))

  # All should have same structure
  expect_true(inherits(fc_link, "tbl_df"))
  expect_true(inherits(fc_expected, "tbl_df"))
  expect_true(inherits(fc_response, "tbl_df"))

  # But estimates should be different
  expect_false(identical(fc_link$.estimate, fc_expected$.estimate))
  expect_false(identical(fc_link$.estimate, fc_response$.estimate))
})

test_that("forecast.ffc_gam() handles summary = FALSE correctly", {
  newdata <- data.frame(
    season = c(1, 2),
    time = c(76, 77)
  )

  fc_raw <- SW(forecast(example_mod,
                        newdata = newdata,
                        summary = FALSE))

  # Should return distributional object
  expect_true(inherits(fc_raw, "distribution"))
  expect_equal(length(fc_raw), nrow(newdata))
})

test_that("forecast.ffc_gam() handles robust vs non-robust summaries", {
  newdata <- data.frame(
    season = c(1, 2),
    time = c(76, 77)
  )

  # Test robust = TRUE (default)
  fc_robust <- SW(forecast(example_mod,
                           newdata = newdata,
                           robust = TRUE))

  # Test robust = FALSE
  fc_nonrobust <- SW(forecast(example_mod,
                              newdata = newdata,
                              robust = FALSE))

  # Both should have same structure
  expect_true(all(c(".estimate", ".error") %in% names(fc_robust)))
  expect_true(all(c(".estimate", ".error") %in% names(fc_nonrobust)))

  # But estimates may be different (median vs mean)
  expect_equal(nrow(fc_robust), nrow(fc_nonrobust))
})

test_that("forecast.ffc_gam() validates type parameter", {
  newdata <- data.frame(
    season = c(1),
    time = c(76)
  )

  # Test invalid type
  expect_error(forecast(example_mod,
                        newdata = newdata,
                        type = "invalid"))
})

test_that("forecast.ffc_gam() handles custom quantiles", {
  newdata <- data.frame(
    season = c(1),
    time = c(76)
  )

  custom_probs <- c(0.1, 0.5, 0.9)
  fc <- SW(forecast(example_mod,
                    newdata = newdata,
                    probs = custom_probs))

  # Check that custom quantile columns exist
  expect_true(all(paste0(".q", 100 * custom_probs) %in% names(fc)))
  expect_equal(ncol(fc), 2 + length(custom_probs))  # .estimate, .error, + quantiles
})

test_that("forecast methods preserve object attributes", {
  functional_coefs <- fts_coefs(example_mod, summary = FALSE, times = 2)

  # Check that fts_ts object has expected attributes
  expect_true(!is.null(attr(functional_coefs, "time_var")))
  expect_true(!is.null(attr(functional_coefs, "summarized")))

  # Forecast should complete without error despite attributes
  fc <- SW(forecast(functional_coefs, h = 1, times = 2))
  expect_true(inherits(fc, "tbl_df"))
})

test_that("forecast handles edge cases appropriately", {
  functional_coefs <- fts_coefs(example_mod, summary = FALSE, times = 2)

  # Test h = 1 (minimum horizon)
  fc_h1 <- SW(forecast(functional_coefs, h = 1, times = 2))
  expect_true(length(unique(fc_h1$.time)) == 1)

  # Test times = 1 (minimum realisations)
  fc_t1 <- SW(forecast(functional_coefs, h = 2, times = 1))
  expect_true(all(fc_t1$.rep == 1))

  # Test larger horizon
  fc_large <- SW(forecast(functional_coefs, h = 5, times = 2))
  expect_true(length(unique(fc_large$.time)) == 5)
})

# Tests for as_fable.ffc_gam function
test_that("as_fable.ffc_gam() creates valid fable objects", {
  skip_if_not_installed("fabletools")
  skip_if_not_installed("distributional")
  skip_if_not_installed("tsibble")

  # Create test newdata
  newdata <- data.frame(
    y = c(10, 15, 20),
    season = c(1, 2, 3),
    time = c(76, 77, 78)
  )

  # Test basic conversion
  fc_fable <- SW(as_fable(example_mod, newdata = newdata))

  # Check fable class structure
  expect_true(inherits(fc_fable, "fbl_ts"))
  expect_true(inherits(fc_fable, "tbl_ts"))
  expect_true(inherits(fc_fable, "tbl_df"))

  # Check required fable columns
  expect_true(all(c(".dist", ".mean", ".model") %in% names(fc_fable)))

  # Check distributional column
  expect_true(inherits(fc_fable$.dist, "distribution"))
  expect_true(is.numeric(fc_fable$.mean))

  # Check model information
  expect_true(all(grepl("FFC_", fc_fable$.model)))

  # Check dimensions
  expect_equal(nrow(fc_fable), nrow(newdata))
})

test_that("as_fable.ffc_gam() works with pre-computed forecasts", {
  skip_if_not_installed("fabletools")
  skip_if_not_installed("distributional")
  skip_if_not_installed("tsibble")

  newdata <- data.frame(
    y = c(10, 15),
    season = c(1, 2),
    time = c(76, 77)
  )

  # Pre-compute forecasts as matrix
  forecasts <- SW(forecast(example_mod, newdata = newdata, summary = FALSE))

  # Convert to fable using pre-computed forecasts
  fc_fable <- SW(as_fable(example_mod, newdata = newdata,
                          forecasts = forecasts))

  expect_true(inherits(fc_fable, "fbl_ts"))
  expect_equal(nrow(fc_fable), nrow(newdata))
  expect_true(inherits(fc_fable$.dist, "distribution"))
})

test_that("as_fable.ffc_gam() auto-detects response variable", {
  skip_if_not_installed("fabletools")
  skip_if_not_installed("distributional")
  skip_if_not_installed("tsibble")

  newdata <- data.frame(
    y = c(10),
    season = c(1),
    time = c(76)
  )

  # Should auto-detect "y" as response
  fc_fable <- SW(as_fable(example_mod, newdata = newdata))

  expect_true(inherits(fc_fable, "fbl_ts"))
  expect_equal(attr(fc_fable, "response"), "y")
})

test_that("as_fable.ffc_gam() handles different forecast models", {
  skip_if_not_installed("fabletools")
  skip_if_not_installed("distributional")
  skip_if_not_installed("tsibble")

  newdata <- data.frame(
    y = c(10),
    season = c(1),
    time = c(76)
  )

  # Test with ETS model
  fc_ets <- SW(as_fable(example_mod, newdata = newdata, model = "ETS"))
  expect_true("FFC_ETS" %in% fc_ets$.model)

  # Test with ARIMA model (default)
  fc_arima <- SW(as_fable(example_mod, newdata = newdata))
  expect_true("FFC_ARIMA" %in% fc_arima$.model)
})

test_that("as_fable.ffc_gam() validates inputs properly", {
  newdata <- data.frame(
    y = c(10),
    season = c(1),
    time = c(76)
  )

  # Test invalid object class
  expect_error(as_fable("not_an_ffc_gam", newdata = newdata))

  # Test missing newdata
  expect_error(as_fable(example_mod, newdata = NULL))

  # Test empty newdata
  expect_error(as_fable(example_mod, newdata = data.frame()))

  # Test missing response in newdata
  bad_newdata <- data.frame(season = c(1), time = c(76))
  expect_error(as_fable(example_mod, newdata = bad_newdata))

  # Test missing time variable in newdata
  bad_newdata2 <- data.frame(y = c(10), season = c(1))
  expect_error(as_fable(example_mod, newdata = bad_newdata2))
})

test_that("as_fable.ffc_gam() handles key variables correctly", {
  skip_if_not_installed("fabletools")
  skip_if_not_installed("distributional")
  skip_if_not_installed("tsibble")

  # Test with grouping variables
  newdata <- data.frame(
    y = c(10, 15, 20, 25),
    season = c(1, 2, 3, 4),
    time = c(76, 77, 78, 79),
    group = c("A", "A", "B", "B")
  )

  fc_fable <- SW(as_fable(example_mod, newdata = newdata))

  # Should be a tsibble with proper key structure
  expect_true(inherits(fc_fable, "tbl_ts"))
  expect_equal(nrow(fc_fable), nrow(newdata))
})

test_that("as_fable.ffc_gam() handles custom key variables", {
  skip_if_not_installed("fabletools")
  skip_if_not_installed("distributional")
  skip_if_not_installed("tsibble")

  newdata <- data.frame(
    y = c(10, 15),
    season = c(1, 2),
    time = c(76, 77),
    custom_key = c("X", "Y")
  )

  # Specify custom key variables
  fc_fable <- SW(as_fable(example_mod, newdata = newdata,
                          key_vars = "custom_key"))

  expect_true(inherits(fc_fable, "fbl_ts"))
  expect_equal(nrow(fc_fable), nrow(newdata))
})

test_that("as_fable.ffc_gam() handles different forecast formats", {
  skip_if_not_installed("fabletools")
  skip_if_not_installed("distributional")
  skip_if_not_installed("tsibble")

  newdata <- data.frame(
    y = c(10, 15),
    season = c(1, 2),
    time = c(76, 77)
  )

  # Test with matrix forecasts
  forecasts_matrix <- matrix(c(10, 12, 14, 16), nrow = 2, ncol = 2)
  fc_matrix <- SW(as_fable(example_mod, newdata = newdata,
                           forecasts = forecasts_matrix))
  expect_true(inherits(fc_matrix, "fbl_ts"))

  # Test with numeric forecasts
  forecasts_numeric <- c(10, 15)
  fc_numeric <- SW(as_fable(example_mod, newdata = newdata,
                            forecasts = forecasts_numeric))
  expect_true(inherits(fc_numeric, "fbl_ts"))
})

test_that("as_fable.ffc_gam() handles binomial cbind() responses", {
  skip_if_not_installed("distributional")

  # Create a binomial model with cbind() response
  binom_data <- data.frame(
    successes = c(8, 6, 4, 5, 9, 7, 9, 8, 9),
    failures = c(2, 4, 6, 5, 1, 3, 3, 1, 5),
    x = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    time = 1:9
  )

  binom_mod <- SW(ffc_gam(
    cbind(successes, failures) ~ fts(time, mean_only = TRUE, time_k = 3, k = 3),
    time = "time",
    data = binom_data,
    family = binomial()
  ))

  # Test auto-detection of cbind() response
  newdata <- data.frame(
    successes = c(7, 8),
    failures = c(3, 2),
    time = c(10, 11)
  )

  fc_fable <- SW(as_fable(binom_mod, newdata = newdata))

  # Should auto-detect cbind response
  expect_true(inherits(fc_fable, "fbl_ts"))
  expected_response <- "cbind(successes,failures)"
  expect_equal(attr(fc_fable, "response"), expected_response)

  # Should handle missing response variables properly
  bad_newdata <- data.frame(
    successes = c(7, 8),
    time = c(7, 8)
  )
  expect_error(as_fable(binom_mod, newdata = bad_newdata))
})

test_that("as_fable.ffc_gam() uses tsibble::new_tsibble correctly", {
  skip_if_not_installed("distributional")

  newdata <- data.frame(
    y = c(10, 15),
    season = c(1, 2),
    time = c(76, 77)
  )

  fc_fable <- SW(as_fable(example_mod, newdata = newdata))

  # Test that the fable object has correct structure
  expect_true(inherits(fc_fable, "fbl_ts"))
  expect_true(inherits(fc_fable, "tbl_ts"))

  # Test that required attributes are set
  expect_equal(attr(fc_fable, "response"), "y")
  expect_equal(attr(fc_fable, "dist"), ".dist")
  expect_equal(attr(fc_fable, "model_cn"), ".model")

  # Test that dimnames are set correctly for distribution
  expect_equal(dimnames(fc_fable$.dist)[[1]], "y")

  # Test distribution column is properly formatted
  expect_true(inherits(fc_fable$.dist, "distribution"))
  expect_true(all(c(".dist", ".mean", ".model") %in% names(fc_fable)))
})
