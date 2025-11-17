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

  # Test ENS ensemble model
  fc_ens <- SW(forecast(functional_coefs,
                        model = "ENS",
                        h = 2,
                        times = 2))
  expect_true(inherits(fc_ens, "tbl_df"))
  expect_true(all(c(".basis", ".sim", ".rep") %in% names(fc_ens)))

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

  # Test basic forecasting with ETS model
  fc <- SW(forecast(example_mod,
                    newdata = newdata,
                    model = "ETS",
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
                         type = "link",
                         model = "ETS"))
  fc_expected <- SW(forecast(example_mod,
                             newdata = newdata,
                             type = "expected",
                             model = "ETS"))
  fc_response <- SW(forecast(example_mod,
                             newdata = newdata,
                             type = "response",
                             model = "ETS"))

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
  forecasts <- SW(forecast(example_mod, newdata = newdata, model = "ETS", summary = FALSE))

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

  # Test with ENS ensemble model
  fc_ens <- SW(as_fable(example_mod, newdata = newdata, model = "ENS"))
  expect_true("FFC_ENS" %in% fc_ens$.model)

  # Test with ARIMA model (default)
  fc_arima <- SW(as_fable(example_mod, newdata = newdata, model = "ARIMA"))
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

test_that("as_fable.ffc_gam() error cases are properly handled", {
  skip_if_not_installed("distributional")

  newdata <- data.frame(
    y = c(10, 15),
    season = c(1, 2),
    time = c(76, 77)
  )

  # Test error when formula is NULL and response not provided
  mod_no_formula <- example_mod
  mod_no_formula$formula <- NULL

  expect_error(
    as_fable(mod_no_formula, newdata = newdata, response = NULL),
    "Cannot auto-detect response variable"
  )

  # Test error when time variable not found in model object
  mod_no_time <- example_mod
  mod_no_time$time_var <- NULL

  expect_error(
    as_fable(mod_no_time, newdata = newdata),
    "Time variable not found in model object"
  )

  # Test error when forecast dimensions don't match newdata
  forecasts_wrong_dim <- matrix(c(10, 12, 14), nrow = 1, ncol = 3)

  expect_error(
    as_fable(example_mod, newdata = newdata, forecasts = forecasts_wrong_dim),
    "Forecast dimensions.*do not match newdata rows"
  )

  # Test time variable class validation with pre-computed forecasts
  # Use pre-computed forecasts to bypass forecast generation issues
  newdata_bad_time <- newdata
  newdata_bad_time$time <- as.complex(newdata_bad_time$time)
  forecasts_good <- matrix(c(10, 12, 14, 16), nrow = 2, ncol = 2)

  expect_error(
    as_fable(example_mod, newdata = newdata_bad_time, forecasts = forecasts_good),
    "must be Date.*POSIXct.*yearquarter.*yearmonth.*numeric.*integer"
  )
})

test_that("as_fable.ffc_gam() handles character time variable conversion", {
  skip_if_not_installed("distributional")

  # Use simple data that won't trigger forecast validation errors
  newdata <- data.frame(
    y = c(10),
    season = c(1),
    time = c("76")  # single character time point
  )

  # Pre-compute forecasts to avoid forecast generation issues
  forecasts_good <- matrix(c(10, 12), nrow = 2, ncol = 1)

  # Should warn about conversion and convert to numeric
  expect_warning(
    fc_fable <- as_fable(example_mod, newdata = newdata, forecasts = forecasts_good),
    "Converting character time variable to numeric"
  )

  expect_true(inherits(fc_fable, "fbl_ts"))
  expect_true(is.numeric(fc_fable$time))
})

# Tests for time interval validation
test_that("forecast.ffc_gam() validates time intervals correctly", {
  # Use growth_data which has regular yearly intervals
  mod_growth <- SW(ffc_gam(
    height_cm ~ s(id, bs = "re") +
      fts(age_yr, k = 4, bs = "cr", time_k = 5),
    time = "age_yr",
    data = growth_data,
    family = Gamma()
  ))

  # Test with correct intervals (should work without warnings)
  newdata_good <- data.frame(
    id = "boy_11",
    age_yr = c(16, 17, 18)
  )

  expect_no_warning({
    fc_good <- forecast(mod_growth, newdata = newdata_good, model = "ETS")
  })
  expect_true(inherits(fc_good, "tbl_df"))
  expect_equal(nrow(fc_good), nrow(newdata_good))

  # Test with mismatched intervals (should warn)
  newdata_warn <- data.frame(
    id = "boy_11",
    age_yr = c(16, 16.5, 17, 17.5)
  )

  expect_warning({
    fc_warn <- forecast(mod_growth, newdata = newdata_warn, model = "ETS")
  }, "differ from training interval")
  expect_true(inherits(fc_warn, "tbl_df"))
  expect_equal(nrow(fc_warn), nrow(newdata_warn))
})

test_that("forecast.ffc_gam() handles overlapping time points", {
  # Use growth_data which has training data from age 3-15
  mod_growth <- SW(ffc_gam(
    height_cm ~ s(id, bs = "re") +
      fts(age_yr, k = 4, bs = "cr", time_k = 5),
    time = "age_yr",
    data = growth_data,
    family = Gamma()
  ))

  # Test with some overlapping time points (should warn and filter)
  newdata_overlap <- data.frame(
    id = "boy_11",
    age_yr = c(14, 15, 16, 17)  # 14, 15 overlap with training data
  )

  expect_warning({
    fc_overlap <- forecast(mod_growth, newdata = newdata_overlap, model = "ETS")
  }, "overlap with training data")

  expect_true(inherits(fc_overlap, "tbl_df"))
  expect_equal(nrow(fc_overlap), 2)  # Should only forecast for 16, 17

  # Test with all overlapping time points (should error)
  newdata_all_overlap <- data.frame(
    id = "boy_11",
    age_yr = c(13, 14, 15)  # All overlap with training data
  )

  expect_error({
    forecast(mod_growth, newdata = newdata_all_overlap, model = "ETS")
  }, "No future time points found")
})

test_that("ENS ensemble works with Quarter time index (tourism example)", {
  # Recreate tourism example from README
  testthat::skip_if_not_installed("lubridate")
  tourism_melb <- tsibble::tourism |>
    dplyr::filter(
      Region == "Melbourne",
      Purpose == "Visiting"
    ) |>
    dplyr::mutate(
      quarter = lubridate::quarter(Quarter),
      time = dplyr::row_number()
    )

  # Split data
  train <- tourism_melb |> dplyr::slice_head(n = 75)
  test <- tourism_melb |> dplyr::slice_tail(n = 5)

  # Fit model
  mod <- SW(ffc_gam(
    Trips ~
      fts(
        time,
        mean_only = TRUE,
        time_k = 50,
        time_m = 1
      ) +
      fts(
        quarter,
        k = 4,
        time_k = 15,
        time_m = 1
      ),
    time = "time",
    data = train,
    family = tw(),
    engine = "gam"
  ))

  # Test ENS ensemble with Quarter time index
  fc_ffc_ens <- SW(as_fable(mod, newdata = test, model = "ENS"))

  # Verify structure
  expect_true(inherits(fc_ffc_ens, "fbl_ts"))
  expect_true("FFC_ENS" %in% fc_ffc_ens$.model)
  expect_true("Quarter" %in% names(fc_ffc_ens))
  expect_equal(nrow(fc_ffc_ens), nrow(test))

  # Test other models still work
  fc_ffc_ets <- SW(as_fable(mod, newdata = test, model = "ETS"))
  fc_ffc_rw <- SW(as_fable(mod, newdata = test, model = "RW"))

  expect_true("FFC_ETS" %in% fc_ffc_ets$.model)
  expect_true("FFC_RW" %in% fc_ffc_rw$.model)
})

test_that("forecast.ffc_gam() handles grouped data time validation", {
  # Create test data with multiple groups
  test_data <- data.frame(
    y = rnorm(30),
    group = rep(c("A", "B", "C"), each = 10),
    time = rep(1:10, 3)
  )

  mod_grouped <- SW(ffc_gam(
    y ~ s(group, bs = "re") + fts(time, k = 4, time_k = 5),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  # Test with correct intervals for all groups
  newdata_good <- data.frame(
    group = c("A", "A", "B", "B"),
    time = c(11, 12, 11, 12)
  )

  expect_no_warning({
    fc_good <- forecast(mod_grouped, newdata = newdata_good, model = "RW")
  })
  expect_true(inherits(fc_good, "tbl_df"))
  expect_equal(nrow(fc_good), nrow(newdata_good))

  # Test with mismatched intervals (should warn)
  newdata_bad <- data.frame(
    group = c("A", "A", "A"),
    time = c(11, 11.5, 12)  # 0.5 intervals vs training 1.0 intervals
  )

  expect_warning({
    fc_bad <- forecast(mod_grouped, newdata = newdata_bad, model = "RW")
  }, "differ from training interval")
})

test_that("adjust_forecast_uncertainty works correctly", {
  # Create mock forecast data that mimics the structure expected by the function
  forecast_df <- data.frame(
    .basis = rep(c("fts_1", "fts_2"), each = 2),
    .realisation = rep(1:2, 2),
    .estimate = distributional::dist_normal(
      mean = c(10, 15, 12, 18),
      sd = rep(0.5, 4)
    ),
    stringsAsFactors = FALSE
  )

  # Create mock standard deviations data
  object_sds <- data.frame(
    .basis = rep(c("fts_1", "fts_2"), each = 2),
    .realisation = rep(1:2, 2),
    .sd = c(1.2, 1.5, 0.8, 1.1),
    stringsAsFactors = FALSE
  )

  times <- 2
  h <- 1

  # Test the internal function
  result <- SW(ffc:::adjust_forecast_uncertainty(
    forecast_df = forecast_df,
    object_sds = object_sds,
    times = times,
    h = h
  ))

  # Check that result has expected structure
  expect_true(inherits(result, "tbl_df"))
  expect_true(all(c(".basis", ".realisation", ".sim", ".rep") %in% names(result)))

  # Check dimensions - should have times * h rows for each basis/realisation combination
  expected_rows <- nrow(forecast_df) * times * h
  expect_equal(nrow(result), expected_rows)

  # Check that .rep variable is properly structured
  expect_true(all(result$.rep %in% 1:times))
  expect_equal(length(unique(result$.rep)), times)

  # Check that .sim values are numeric
  expect_true(is.numeric(result$.sim))
  expect_true(all(is.finite(result$.sim)))

  # Check that all original basis and realisation combinations are preserved
  original_combos <- paste(forecast_df$.basis, forecast_df$.realisation)
  result_combos <- unique(paste(result$.basis, result$.realisation))
  expect_setequal(result_combos, original_combos)

  # Check that standard deviations were properly joined
  expect_true(".sd" %in% names(result))
  expect_true(all(is.finite(result$.sd)))
})

test_that("adjust_forecast_uncertainty handles missing .sd data correctly", {
  # Create forecast data
  forecast_df <- data.frame(
    .basis = c("fts_1", "fts_2"),
    .realisation = c(1, 1),
    .estimate = distributional::dist_normal(mean = c(10, 15), sd = c(0.5, 0.7)),
    stringsAsFactors = FALSE
  )

  # Create object_sds with only partial matches
  object_sds <- data.frame(
    .basis = "fts_1",
    .realisation = 1,
    .sd = 1.2,
    stringsAsFactors = FALSE
  )

  times <- 2
  h <- 1

  # Should handle the missing .sd gracefully through the left_join
  result <- SW(ffc:::adjust_forecast_uncertainty(
    forecast_df = forecast_df,
    object_sds = object_sds,
    times = times,
    h = h
  ))

  expect_true(inherits(result, "tbl_df"))
  expect_equal(nrow(result), nrow(forecast_df) * times * h)

  # Check that we have some NA values where .sd wasn't available
  expect_true(any(is.na(result$.sd)))
  expect_true(any(!is.na(result$.sd)))
})

test_that("forecast.ffc_gam() works with models without fts() terms", {
  # Create a model without any fts() terms - just a regular GAM
  test_data <- data.frame(
    y = rnorm(50, mean = 10, sd = 2),
    x1 = runif(50, 0, 10),
    x2 = rnorm(50, 0, 1),
    time = 1:50
  )

  # Fit a model without fts() terms
  mod_no_fts <- SW(ffc_gam(
    y ~ s(x1, k = 5) + s(x2, k = 5),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  # Check that gam_init is empty (no fts terms)
  expect_equal(length(mod_no_fts$gam_init), 0)

  # Create newdata for forecasting
  newdata <- data.frame(
    x1 = runif(5, 0, 10),
    x2 = rnorm(5, 0, 1),
    time = 51:55
  )

  # Test that forecast works and returns predictions
  fc_no_fts <- SW(forecast(mod_no_fts, newdata = newdata))

  # Check that result has expected structure
  expect_true(inherits(fc_no_fts, "tbl_df"))
  expect_true(all(c(".estimate", ".error") %in% names(fc_no_fts)))
  expect_equal(nrow(fc_no_fts), nrow(newdata))

  # Check that predictions are numeric and finite
  expect_true(is.numeric(fc_no_fts$.estimate))
  expect_true(all(is.finite(fc_no_fts$.estimate)))
  expect_true(is.numeric(fc_no_fts$.error))
  expect_true(all(fc_no_fts$.error >= 0))

  # Test with different types
  fc_link <- SW(forecast(mod_no_fts, newdata = newdata, type = "link"))
  fc_response <- SW(forecast(mod_no_fts, newdata = newdata, type = "response"))

  expect_true(inherits(fc_link, "tbl_df"))
  expect_true(inherits(fc_response, "tbl_df"))

  # Check that both produce numeric results (may differ due to link function)
  expect_true(is.numeric(fc_link$.estimate))
  expect_true(is.numeric(fc_response$.estimate))
  expect_true(all(is.finite(fc_link$.estimate)))
  expect_true(all(is.finite(fc_response$.estimate)))
})

test_that("forecast.ffc_gam() works with mixed models (some fts, some regular smooths)", {
  # Create a model with both fts() and regular smooth terms
  test_data <- data.frame(
    y = rnorm(60, mean = 15, sd = 3),
    x1 = runif(60, 0, 10),
    x2 = rnorm(60, 0, 2),
    time = 1:60
  )

  # Fit a model with mixed terms
  mod_mixed <- SW(ffc_gam(
    y ~ s(x1, k = 5) + fts(x2, k = 4, time_k = 5),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  # Check that the model has fts_smooths
  expect_false(is.null(mod_mixed$fts_smooths))
  expect_true(length(mod_mixed$fts_smooths) > 0)
  expect_true(length(mod_mixed$gam_init) > 0)

  # Create newdata for forecasting
  newdata <- data.frame(
    x1 = runif(5, 0, 10),
    x2 = rnorm(5, 0, 2),
    time = 61:65
  )

  # This should work since we have fts terms to forecast
  fc_mixed <- SW(forecast(mod_mixed, newdata = newdata))

  expect_true(inherits(fc_mixed, "tbl_df"))
  expect_true(all(c(".estimate", ".error") %in% names(fc_mixed)))
  expect_equal(nrow(fc_mixed), nrow(newdata))

  # Test that coefficients can be extracted from mixed model
  coefs <- fts_coefs(mod_mixed, summary = TRUE)
  expect_true(inherits(coefs, "fts_ts"))

  # Coefficients should only be for the fts terms, not regular smooths
  expect_true(all(grepl("fts", coefs$.basis)))
})

test_that("adjust_forecast_uncertainty preserves group structure", {
  # Create forecast data with multiple groups
  forecast_df <- data.frame(
    .basis = c("fts_1", "fts_1", "fts_2", "fts_2"),
    .realisation = c(1, 2, 1, 2),
    .estimate = distributional::dist_normal(
      mean = c(2, 4, 6, 8),
      sd = rep(0.3, 4)
    ),
    stringsAsFactors = FALSE
  )

  object_sds <- data.frame(
    .basis = c("fts_1", "fts_1", "fts_2", "fts_2"),
    .realisation = c(1, 2, 1, 2),
    .sd = c(0.8, 1.0, 0.6, 1.2),
    stringsAsFactors = FALSE
  )

  times <- 2
  h <- 1

  result <- SW(ffc:::adjust_forecast_uncertainty(
    forecast_df = forecast_df,
    object_sds = object_sds,
    times = times,
    h = h
  ))

  # Check that each basis/realisation group has the correct number of rows
  group_counts <- result |>
    dplyr::group_by(.basis, .realisation) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  expect_true(all(group_counts$n == times * h))

  # Check that .rep is correctly structured within each group
  rep_structure <- result |>
    dplyr::group_by(.basis, .realisation) |>
    dplyr::summarise(
      min_rep = min(.rep),
      max_rep = max(.rep),
      unique_reps = length(unique(.rep)),
      .groups = "drop"
    )

  expect_true(all(rep_structure$min_rep == 1))
  expect_true(all(rep_structure$max_rep == times))
  expect_true(all(rep_structure$unique_reps == times))
})
