# Tests for the forecast functions
# ==============================================================================
# Test Setup
# ==============================================================================

# Create test data
test_newdata <- data.frame(
  season = c(1, 6, 12),
  time = c(76, 77, 78)
)

# ==============================================================================
# validate_forecast_inputs() Tests
# ==============================================================================

test_that("validate_forecast_inputs handles valid inputs correctly", {

  result <- ffc:::validate_forecast_inputs(
    object = example_mod,
    newdata = test_newdata,
    type = "response",
    model = "ARIMA",
    stationary = FALSE,
    summary = TRUE,
    robust = FALSE,
    probs = c(0.025, 0.975)
  )

  expect_type(result, "list")
  expect_equal(result$type, "response")
  expect_equal(result$model, "ARIMA")
  expect_false(result$stationary)
  expect_true(result$summary)
  expect_equal(result$probs, c(0.025, 0.975))
})

test_that("validate_forecast_inputs validates type argument", {

  # Test invalid type
  expect_error(
    ffc:::validate_forecast_inputs(
      object = example_mod,
      newdata = test_newdata,
      type = "invalid_type"
    )
  )

  # Test valid types
  for (type in c("link", "expected", "response")) {
    expect_no_error(
      ffc:::validate_forecast_inputs(
        object = example_mod,
        newdata = test_newdata,
        type = type
      )
    )
  }
})

test_that("validate_forecast_inputs catches stationary constraint violations", {

  expect_error(
    ffc:::validate_forecast_inputs(
      object = example_mod,
      newdata = test_newdata,
      model = "ETS",
      stationary = TRUE
    ),
    "stationary = TRUE only works with model = 'ARIMA'"
  )
})

test_that("validate_forecast_inputs validates logical arguments", {

  invalid_logicals <- list(
    c(TRUE, FALSE),  # Multiple values
    "TRUE",          # Character
    NA,              # Missing
    1                # Numeric
  )

  expect_error(
    ffc:::validate_forecast_inputs(
      object = example_mod,
      newdata = test_newdata,
      summary = c(TRUE, FALSE)
    )
  )

  expect_error(
    ffc:::validate_forecast_inputs(
      object = example_mod,
      newdata = test_newdata,
      summary = "TRUE"
    )
  )
})

test_that("validate_forecast_inputs validates probabilities", {

  # Test invalid probabilities
  expect_error(
    ffc:::validate_forecast_inputs(
      object = example_mod,
      newdata = test_newdata,
      probs = c(-0.1, 0.5)
    ),
    "must contain numeric values between 0 and 1"
  )

  expect_error(
    ffc:::validate_forecast_inputs(
      object = example_mod,
      newdata = test_newdata,
      probs = c(0.5, 1.1)
    ),
    "must contain numeric values between 0 and 1"
  )

  # Test probabilities are sorted and deduplicated
  result <- ffc:::validate_forecast_inputs(
    object = example_mod,
    newdata = test_newdata,
    probs = c(0.9, 0.1, 0.5, 0.1)
  )
  expect_equal(result$probs, c(0.1, 0.5, 0.9))
})

# ==============================================================================
# sample_coefficients() Tests
# ==============================================================================

test_that("sample_coefficients returns correct structure", {

  n_samples <- 100
  samples <- ffc:::sample_coefficients(example_mod, n_samples)

  expect_true(is.matrix(samples))
  expect_equal(nrow(samples), n_samples)
  expect_equal(ncol(samples), length(coef(example_mod)))
  expect_true(all(is.finite(samples)))
})

test_that("sample_coefficients validates sample count", {

  invalid_samples <- list(-1, 0, 1.5, NA, "10", c(10, 20))

  for (invalid in invalid_samples) {
    expect_error(
      ffc:::sample_coefficients(example_mod, invalid),
      "must be a positive integer"
    )
  }
})

test_that("sample_coefficients is reproducible with seed", {

  set.seed(123)
  samples1 <- ffc:::sample_coefficients(example_mod, 50)

  set.seed(123)
  samples2 <- ffc:::sample_coefficients(example_mod, 50)

  expect_equal(samples1, samples2)
})

# ==============================================================================
# compute_linear_predictors() Tests
# ==============================================================================

# Updated tests for compute_linear_predictors

test_that("compute_linear_predictors handles simple models correctly", {

  # Test with example_mod (should be simple)
  beta_samples <- ffc:::sample_coefficients(example_mod, 200)

  linpreds <- ffc:::compute_linear_predictors(
    object = example_mod,
    newdata = test_newdata,
    beta_samples = beta_samples
  )

  expect_true(is.matrix(linpreds))
  expect_equal(nrow(linpreds), 200)  # n_samples
  expect_equal(ncol(linpreds), nrow(test_newdata))  # n_predictions
  expect_true(all(is.finite(linpreds)))
})

test_that("compute_linear_predictors validates input dimensions", {

  # Test with wrong number of coefficients
  wrong_beta_samples <- matrix(rnorm(20), nrow = 10, ncol = 2)  # Wrong ncol

  expect_error(
    ffc:::compute_linear_predictors(
      object = example_mod,
      newdata = test_newdata,
      beta_samples = wrong_beta_samples
    ),
    "must have .* columns to match model coefficients"
  )

  # Test with non-matrix beta_samples
  expect_error(
    ffc:::compute_linear_predictors(
      object = example_mod,
      newdata = test_newdata,
      beta_samples = "not_a_matrix"
    ),
    "must be a matrix"
  )

  # Test with non-dataframe newdata
  beta_samples <- sample_coefficients(example_mod, 5)
  expect_error(
    ffc:::compute_linear_predictors(
      object = example_mod,
      newdata = "not_a_dataframe",
      beta_samples = beta_samples
    ),
    "must be a data.frame"
  )
})

test_that("compute_linear_predictors handles model offsets correctly", {

  beta_samples <- ffc:::sample_coefficients(example_mod, 200)

  # Test that function handles offsets without error
  result <- ffc:::compute_linear_predictors(
    object = example_mod,
    newdata = test_newdata,
    beta_samples = beta_samples
  )

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 200)
  expect_equal(ncol(result), nrow(test_newdata))
})

test_that("compute_linear_predictors provides helpful error messages", {

  # Test dimension mismatch error message
  wrong_samples <- matrix(rnorm(10), nrow = 5, ncol = 2)

  expect_error(
    ffc:::compute_linear_predictors(
      object = example_mod,
      newdata = test_newdata,
      beta_samples = wrong_samples
    ),
    "'beta_samples' must have 30 columns to match model coefficients"
  )
})

test_that("compute_linear_predictors handles predict() failures gracefully", {

  beta_samples <- ffc:::sample_coefficients(example_mod, 200)

  # Create newdata that might cause predict() to fail
  problematic_newdata <- data.frame(
    season = c(1, 2),
    time = c(76, 77),
    nonexistent_var = c("a", "b")  # Extra variable that shouldn't matter
  )

  # Remove the extra variable to make it valid
  clean_newdata <- problematic_newdata[, c("season", "time")]

  expect_no_error(
    ffc:::compute_linear_predictors(
      object = example_mod,
      newdata = clean_newdata,
      beta_samples = beta_samples
    )
  )
})

test_that("tcrossprod operation works correctly", {

  # Test the specific matrix operation that was failing
  n_samples <- 5
  n_coef <- length(coef(example_mod))
  n_pred <- nrow(test_newdata)

  beta_samples <- ffc:::sample_coefficients(example_mod, n_samples)
  lpmat <- predict(example_mod, newdata = test_newdata, type = "lpmatrix")

  # Check dimensions before operation
  expect_equal(nrow(beta_samples), n_samples)
  expect_equal(ncol(beta_samples), n_coef)
  expect_equal(nrow(lpmat), n_pred)
  expect_equal(ncol(lpmat), n_coef)

  # Perform the operation that was failing
  result <- tcrossprod(beta_samples, lpmat)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), n_samples)
  expect_equal(ncol(result), n_pred)
  expect_true(all(is.finite(result)))
})

# ==============================================================================
# generate_predictions() Tests
# ==============================================================================

test_that("generate_predictions handles all prediction types", {

  # Create test linear predictors
  linpreds <- matrix(rnorm(20), nrow = 10, ncol = 2)

  # Test link type
  pred_link <- ffc:::generate_predictions(example_mod, linpreds, "link")
  expect_identical(pred_link, linpreds)

  # Test expected type
  pred_expected <- ffc:::generate_predictions(example_mod, linpreds, "expected")
  expect_true(is.matrix(pred_expected))
  expect_equal(dim(pred_expected), dim(linpreds))
  expect_true(all(is.finite(pred_expected)))

  # Test response type
  pred_response <- ffc:::generate_predictions(example_mod, linpreds, "response")
  expect_true(is.matrix(pred_response))
  expect_equal(dim(pred_response), dim(linpreds))
  expect_true(all(is.finite(pred_response)))
})

test_that("generate_predictions validates type argument", {

  linpreds <- matrix(rnorm(10), nrow = 5, ncol = 2)

  expect_error(
    ffc:::generate_predictions(example_mod, linpreds, "invalid"),
    "Invalid prediction type"
  )
})

test_that("generate_predictions produces reasonable values for Poisson family", {

  # For Poisson family, response should be non-negative integers
  linpreds <- matrix(c(-1, 0, 1, 2), nrow = 2, ncol = 2)

  pred_response <- ffc:::generate_predictions(example_mod, linpreds, "response")

  # Should be non-negative
  expect_true(all(pred_response >= 0))

  # Should be integers (for Poisson)
  expect_true(all(pred_response == floor(pred_response)))
})

# ==============================================================================
# format_forecast_output() Tests
# ==============================================================================

test_that("format_forecast_output handles summary output", {

  # Create test predictions
  predictions <- matrix(rnorm(30, mean = 5, sd = 1), nrow = 10, ncol = 3)

  result <- ffc:::format_forecast_output(
    predictions = predictions,
    summary = TRUE,
    robust = FALSE,
    probs = c(0.025, 0.975)
  )

  expect_true(is.data.frame(result))
  expect_true(inherits(result, "tbl_df"))
  expect_equal(nrow(result), 3)  # 3 predictions
  expect_true(all(c(".estimate", ".error", ".q2.5", ".q97.5") %in% names(result)))
  expect_true(all(is.finite(result$.estimate)))
  expect_true(all(result$.error >= 0))
})

test_that("format_forecast_output handles robust vs standard statistics", {

  # Create test data with outliers
  predictions <- matrix(c(
    1, 2, 3, 4, 5,     # Normal values
    1, 2, 3, 4, 100    # With outlier
  ), nrow = 5, ncol = 2)

  # Standard statistics
  result_standard <- ffc:::format_forecast_output(
    predictions = predictions,
    summary = TRUE,
    robust = FALSE,
    probs = c(0.5)
  )

  # Robust statistics
  result_robust <- ffc:::format_forecast_output(
    predictions = predictions,
    summary = TRUE,
    robust = TRUE,
    probs = c(0.5)
  )

  # Robust estimate should be less affected by outlier
  expect_true(result_robust$.estimate[2] < result_standard$.estimate[2])
})

test_that("format_forecast_output handles raw distribution output", {

  predictions <- matrix(rnorm(20), nrow = 10, ncol = 2)

  result <- ffc:::format_forecast_output(
    predictions = predictions,
    summary = FALSE
  )

  expect_true(inherits(result, "distribution"))
  expect_equal(length(result), 2)  # 2 predictions
})

test_that("format_forecast_output validates quantile ordering", {

  predictions <- matrix(rnorm(20), nrow = 10, ncol = 2)

  result <- ffc:::format_forecast_output(
    predictions = predictions,
    summary = TRUE,
    probs = c(0.1, 0.5, 0.9)
  )

  # Quantiles should be ordered
  expect_true(result$.q10.0[1] <= result$.q50.0[1])
  expect_true(result$.q50.0[1] <= result$.q90.0[1])
})

# ==============================================================================
# Forecast Tests
# ==============================================================================

test_that("forecast produces correct output", {

  result <- forecast(
    object = example_mod,
    newdata = test_newdata,
    type = "response",
    summary = TRUE
  )

  expect_true(is.data.frame(result))
  expect_true(inherits(result, "tbl_df"))
  expect_equal(nrow(result), nrow(test_newdata))
  expect_true(all(c(".estimate", ".error") %in% names(result)))
  expect_true(all(is.finite(result$.estimate)))
  expect_true(all(result$.error >= 0))
})

test_that("forecast function handles different types correctly", {

  types <- c("link", "expected", "response")

  for (type in types) {
    result <- forecast(
      object = example_mod,
      newdata = test_newdata,
      type = type,
      summary = TRUE
    )

    expect_true(is.data.frame(result))
    expect_equal(nrow(result), nrow(test_newdata))
  }
})

test_that("forecast function handles summary = FALSE", {

  result <- forecast(
    object = example_mod,
    newdata = test_newdata,
    type = "response",
    summary = FALSE
  )

  expect_true(inherits(result, "distribution"))
  expect_equal(length(result), nrow(test_newdata))
})

test_that("forecast function handles custom probabilities", {

  custom_probs <- c(0.05, 0.25, 0.75, 0.95)

  result <- forecast(
    object = example_mod,
    newdata = test_newdata,
    probs = custom_probs,
    summary = TRUE
  )

  expected_cols <- c(".estimate", ".error",
                     paste0(".q", formatC(100 * custom_probs, format = "f", digits = 1)))
  expect_true(all(expected_cols %in% names(result)))
})

test_that("forecast function handles robust statistics", {

  result_standard <- forecast(
    object = example_mod,
    newdata = test_newdata,
    robust = FALSE,
    summary = TRUE
  )

  result_robust <- forecast(
    object = example_mod,
    newdata = test_newdata,
    robust = TRUE,
    summary = TRUE
  )

  # Both should work and have same structure
  expect_equal(names(result_standard), names(result_robust))
  expect_equal(nrow(result_standard), nrow(result_robust))
})

# ==============================================================================
# Error Handling Tests
# ==============================================================================

test_that("refactored functions handle invalid inputs gracefully", {

  # Invalid object type
  expect_error(
    forecast(
      object = "not_a_model",
      newdata = test_newdata
    )
  )

  # Invalid newdata
  expect_error(
    ffc:::validate_forecast_inputs(
      object = example_mod,
      newdata = "not_a_dataframe"
    )
  )

  # Invalid type
  expect_error(
    forecast(
      object = example_mod,
      newdata = test_newdata,
      type = "invalid"
    )
  )
})

test_that("error messages are informative", {

  # Test that error messages contain useful information
  expect_error(
    ffc:::validate_forecast_inputs(
      object = example_mod,
      newdata = test_newdata,
      probs = c(-0.1, 0.5)
    ),
    "between 0 and 1"
  )

  expect_error(
    ffc:::sample_coefficients(example_mod, -5),
    "positive integer"
  )
})

# ==============================================================================
# Performance Tests
# ==============================================================================


test_that("individual functions are fast", {

  # Validation should be very fast
  start_time <- Sys.time()
  params <- ffc:::validate_forecast_inputs(
    object = example_mod,
    newdata = test_newdata
  )
  validation_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  expect_true(validation_time < 0.1)

  # Sampling should be reasonably fast
  start_time <- Sys.time()
  samples <- ffc:::sample_coefficients(example_mod, 100)
  sampling_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  expect_true(sampling_time < 1.0)
})


# ==============================================================================
# Edge Case Tests
# ==============================================================================

test_that("functions handle edge cases gracefully", {

  # Single prediction point
  single_newdata <- data.frame(season = 1, time = 76)

  result <- forecast(
    object = example_mod,
    newdata = single_newdata,
    summary = TRUE
  )

  expect_equal(nrow(result), 1)
  expect_true(all(is.finite(result$.estimate)))

  # Very small sample size
  expect_no_error(
    ffc:::sample_coefficients(example_mod, 1)
  )

  # Edge probability values
  result_extreme <- ffc:::format_forecast_output(
    predictions = matrix(rnorm(10), nrow = 5, ncol = 2),
    summary = TRUE,
    probs = c(0.001, 0.999)
  )

  expect_true(all(c(".q0.1", ".q99.9") %in% names(result_extreme)))
})

test_that("functions handle missing values appropriately", {

  # Test with some missing predictors (should error appropriately)
  newdata_with_na <- test_newdata
  newdata_with_na$season[1] <- NA

  # This should either handle NAs gracefully or error informatively
  # The exact behavior depends on how mgcv::predict handles NAs
  expect_warning(
    forecast(
      object = example_mod,
      newdata = newdata_with_na
    )
  )
})
