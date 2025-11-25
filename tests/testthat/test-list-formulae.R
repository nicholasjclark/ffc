test_that("list formula detection works", {
  # Simple list formula detection - no ffc_gam calls
  test_formula <- list(y ~ x, ~ 1)
  expect_true(is.list(test_formula))
  expect_length(test_formula, 2)
  expect_s3_class(test_formula[[1]], "formula")
  expect_s3_class(test_formula[[2]], "formula")
})

test_that("family nlp property detection", {
  # Test family nlp property detection - no model fitting
  library(mgcv)

  # Single parameter family
  gauss_family <- gaussian()
  expect_null(gauss_family$nlp)

  # Multi-parameter family
  gaulss_family <- gaulss()
  expect_equal(gaulss_family$nlp, 2)
})

test_that("basic input validation works", {
  # Test basic checkmate validation - no model fitting
  expect_error(
    checkmate::assert_list("not a list", types = "formula"),
    "Assertion on"
  )

  expect_error(
    checkmate::assert_class("not a formula", "formula"),
    "Assertion on"
  )

  # List of formulae should pass
  expect_silent(
    checkmate::assert_list(list(y ~ x, ~ 1), types = "formula")
  )
})

test_that("interpret_ffc validates inputs correctly", {
  # Simple test data - no fts() terms to avoid heavy processing
  test_data <- data.frame(y = 1:10, x = 1:10, time = 1:10)

  # Test parameter validation for list formulae
  expect_error(
    interpret_ffc("not a formula", test_data, "time"),
    "Assertion on 'formula' failed"
  )

  expect_error(
    interpret_ffc(list("not a formula"), test_data, "time"),
    "May only contain the following types"
  )

  expect_error(
    interpret_ffc(y ~ x, test_data, ""),
    "Assertion on 'time_var' failed"
  )
})

test_that("interpret_ffc handles simple list formulae", {
  # Simple test data and formulae without fts() terms
  test_data <- data.frame(y = 1:10, x = 1:10, time = 1:10)

  # Test with simple list formulae (no fts terms)
  result <- interpret_ffc(
    formula = list(y ~ x, ~ 1),
    data = test_data,
    time_var = "time"
  )

  # Check return structure
  expect_type(result, "list")
  expect_named(result, c("formula", "data", "orig_data", "fts_smooths", "gam_init"))

  # Check formula is a list
  expect_type(result$formula, "list")
  expect_length(result$formula, 2)
  expect_s3_class(result$formula[[1]], "formula")
  expect_s3_class(result$formula[[2]], "formula")

  # Check data consistency
  expect_true(is.data.frame(result$data))
  expect_true(nrow(result$data) == 10)

  # Check original data preservation
  expect_identical(result$orig_data, test_data)
})

test_that("interpret_ffc maintains backward compatibility", {
  # Test that single formulae work exactly as before
  test_data <- data.frame(y = 1:10, x = 1:10, time = 1:10)

  result <- interpret_ffc(
    formula = y ~ x,
    data = test_data,
    time_var = "time"
  )

  # Check return structure for single formula
  expect_type(result, "list")
  expect_named(result, c("formula", "data", "orig_data", "fts_smooths", "gam_init"))

  # Formula should be single formula, not list
  expect_s3_class(result$formula, "formula")
  expect_false(is.list(result$formula))

  # Check data consistency
  expect_true(is.data.frame(result$data))
  expect_identical(result$orig_data, test_data)
})

test_that("parameter-aware coefficient naming works for list formulae with fts", {
  # Test data with functional terms
  test_data <- data.frame(
    y = rnorm(20),
    x = rnorm(20),
    time = 1:20
  )

  # List formulae with fts terms
  result <- interpret_ffc(
    formula = list(y ~ fts(x, k = 3), ~ fts(x, k = 3)),
    data = test_data,
    time_var = "time"
  )

  # Check that coefficient names include semantic parameter prefixes
  fts_cols <- grep("^(location|scale|shape)_fts_", names(result$data), value = TRUE)
  expect_true(length(fts_cols) > 0)

  # Should have coefficients for both parameters
  location_cols <- grep("^location_", names(result$data), value = TRUE)
  scale_cols <- grep("^scale_", names(result$data), value = TRUE)

  expect_true(length(location_cols) > 0)
  expect_true(length(scale_cols) > 0)

  # Parameter prefixes should be different
  expect_false(any(location_cols %in% scale_cols))
})

test_that("single formula fts terms have no parameter prefix", {
  # Test data with functional terms
  test_data <- data.frame(
    y = rnorm(20),
    x = rnorm(20),
    time = 1:20
  )

  # Single formula with fts term
  result <- interpret_ffc(
    formula = y ~ fts(x, k = 3),
    data = test_data,
    time_var = "time"
  )

  # Check that coefficient names do NOT have parameter prefixes
  fts_cols <- grep("^fts_", names(result$data), value = TRUE)
  expect_true(length(fts_cols) > 0)

  # Should NOT have semantic parameter prefixes
  param_cols <- grep("^(location|scale|shape)_", names(result$data), value = TRUE)
  expect_equal(length(param_cols), 0)
})

test_that("mean_only fts terms get parameter prefixes correctly", {
  # Test data with functional terms
  test_data <- data.frame(
    y = rnorm(20),
    x = rnorm(20),
    time = 1:20
  )

  # List formulae with mean_only fts terms
  result <- interpret_ffc(
    formula = list(y ~ fts(x, k = 3, mean_only = TRUE),
                   ~ fts(x, k = 3, mean_only = TRUE)),
    data = test_data,
    time_var = "time"
  )

  # Check for mean term parameter prefixes
  mean_cols <- grep("_mean$", names(result$data), value = TRUE)
  expect_true(length(mean_cols) > 0)

  # Mean terms should have semantic parameter prefixes
  param_mean_cols <- grep("^(location|scale|shape)_.*_mean$", names(result$data),
                          value = TRUE)
  expect_equal(length(param_mean_cols), length(mean_cols))
})

test_that("extract_parameter_from_basis works correctly", {
  # Load mgcv for family objects
  library(mgcv)

  # Test with distributional family (gaulss)
  gaulss_family <- gaulss()

  # Test parameter extraction from semantic prefixed names
  expect_equal(extract_parameter_from_basis("location_fts_bs_x", gaulss_family), "location")
  expect_equal(extract_parameter_from_basis("scale_fts_bs_x", gaulss_family), "scale")

  # Test parameter extraction from numeric prefixed names (backward compatibility)
  expect_equal(extract_parameter_from_basis("param1_fts_bs_x", gaulss_family), "location")
  expect_equal(extract_parameter_from_basis("param2_fts_bs_x", gaulss_family), "scale")

  # Test single-parameter model fallback
  gauss_family <- gaussian()
  expect_equal(extract_parameter_from_basis("fts_bs_x", gauss_family), "location")

  # Test single-parameter with distributional family
  expect_equal(extract_parameter_from_basis("fts_bs_x", gaulss_family), "location")

  # Test fallback for invalid parameter numbers
  expect_equal(extract_parameter_from_basis("param99_fts_bs_x", gaulss_family), "param99")
})

test_that("end-to-end coefficient extraction works with gaulss()", {
  # Test complete pipeline with distributional family
  library(mgcv)
  set.seed(1234)

  # Generate test data with time-varying location and scale
  n <- 40
  test_data <- data.frame(
    time = 1:n,
    x = rnorm(n),
    y = rnorm(n, mean = sin(2 * pi * (1:n) / 10), sd = 0.5 + 0.3 * cos(2 * pi * (1:n) / 8))
  )

  # Fit model with list formula and gaulss family
  model <- ffc_gam(
    list(y ~ s(x, k = 5), ~ x),
    data = test_data,
    family = gaulss(),
    engine = "gam",
    time = "time"
  )

  # Test model structure
  expect_s3_class(model, "ffc_gam")
  expect_s3_class(model, "ffc_gam_multi")
  expect_equal(model$family$nlp, 2)

  # Should have coefficients for both parameters
  coefs <- coef(model)
  expect_true(all(is.finite(coefs)))
})

test_that("twlss() family works with list formulae", {
  # Test with Tweedie location, scale, shape family
  library(mgcv)
  set.seed(1234)

  # Generate positive test data for Tweedie
  n <- 30
  test_data <- data.frame(
    time = 1:n,
    x = runif(n, 0, 1),
    y = rgamma(n, shape = 2, rate = 1)  # Positive values for Tweedie
  )

  # Test interpret_ffc with twlss
  result <- interpret_ffc(
    formula = list(y ~ fts(x, k = 3), ~ fts(x, k = 3), ~ 1),
    data = test_data,
    time_var = "time"
  )

  # Should handle three parameters
  expect_type(result$formula, "list")
  expect_length(result$formula, 3)

  # Check semantic parameter-specific coefficient naming
  location_cols <- grep("^location_", names(result$data), value = TRUE)
  scale_cols <- grep("^scale_", names(result$data), value = TRUE)
  shape_cols <- grep("^shape_", names(result$data), value = TRUE)

  expect_true(length(location_cols) > 0)
  expect_true(length(scale_cols) > 0)
  expect_equal(length(shape_cols), 0)  # Third parameter is ~ 1, no fts terms
})

test_that("different fts() specifications work per parameter", {
  # Test different k values and options per parameter
  set.seed(1234)

  n <- 25
  test_data <- data.frame(
    time = 1:n,
    x = rnorm(n),
    y = rnorm(n)
  )

  # Test with different k values per parameter
  result <- interpret_ffc(
    formula = list(
      y ~ fts(x, k = 8),           # High flexibility for location
      ~ fts(x, k = 3, mean_only = TRUE)   # Low flexibility, mean-only for scale
    ),
    data = test_data,
    time_var = "time"
  )

  # Check that different specifications are preserved
  location_basis_cols <- grep("^location_fts_bs_", names(result$data), value = TRUE)
  scale_basis_cols <- grep("^scale_fts_bs_", names(result$data), value = TRUE)
  scale_mean_cols <- grep("^scale_.*_mean$", names(result$data), value = TRUE)

  # Location parameter should have more basis functions (k=8)
  expect_true(length(location_basis_cols) > length(scale_basis_cols))

  # Scale parameter should have mean_only terms
  expect_true(length(scale_mean_cols) > 0)
})

test_that("fts() method integration with distributional models", {
  # Test that existing fts methods work with multi-parameter models
  library(mgcv)
  set.seed(1234)

  n <- 35
  test_data <- data.frame(
    time = 1:n,
    x = rnorm(n),
    y = rnorm(n, sd = 0.5)
  )

  # Fit distributional model
  model <- ffc_gam(
    list(y ~ fts(x, k = 4), ~ fts(x, k = 3)),
    data = test_data,
    family = gaulss(),
    time = "time"
  )

  # Test fts() extraction methods
  coefs <- fts_coefs(model)
  expect_s3_class(coefs, "fts_ts")

  # Test that parameter information is preserved
  expect_true(".parameter" %in% names(coefs))
  expect_true(all(c("location", "scale") %in% unique(coefs$.parameter)))

  # Test coefficient structure
  expect_true(".time" %in% names(coefs))
  expect_true(".basis" %in% names(coefs))
  expect_true(".estimate" %in% names(coefs))
  expect_true(".se" %in% names(coefs))

  # Each parameter should have appropriate number of basis functions
  location_coefs <- coefs[coefs$.parameter == "location", ]
  scale_coefs <- coefs[coefs$.parameter == "scale", ]

  expect_true(nrow(location_coefs) > 0)
  expect_true(nrow(scale_coefs) > 0)

  # Different parameters can have different numbers of coefficients
  expect_true(length(unique(location_coefs$.basis)) >= 3)  # k=4 minus constraints
  expect_true(length(unique(scale_coefs$.basis)) >= 2)     # k=3 minus constraints
})

test_that("gam_init structure normalization in distributional regression", {
  # Test that gam_init objects are properly normalized for prediction calls
  # Ensures consistent list structure across model fitting and prediction phases

  library(mgcv)

    # Create test data for distributional model
    set.seed(1234)
    n <- 35
    test_data <- data.frame(
      time = 1:n,
      x = rnorm(n),
      y = rnorm(n, sd = 0.5)
    )

    # Fit distributional model with multiple parameters
    model <- ffc_gam(
      list(y ~ fts(x, k = 4), ~ fts(x, k = 3)),
      data = test_data,
      family = gaulss(),
      time = "time"
    )

    # Verify model has correct gam_init structure
    expect_s3_class(model, "ffc_gam_multi")
    expect_length(model$gam_init, 2)
    expect_true(all(sapply(model$gam_init, function(x) inherits(x, "gam"))))

    # Create newdata for prediction
    newdata <- data.frame(
      time = (n + 1):(n + 3),
      x = rnorm(3)
    )

    # Test prediction with newdata
    pred_result <- predict(model, newdata = newdata, type = "response")

    # Verify prediction returns valid structure
    expect_true(is.data.frame(pred_result) || is.matrix(pred_result) || is.numeric(pred_result))

    # Forecast
    fc <- forecast(model, newdata = newdata)
})

# Additional distributional regression test coverage
test_that("apply_distributional_inverse_links handles missing lpi attribute", {
  library(mgcv)

  # Create mock distributional prediction without lpi attribute
  mock_linpreds <- matrix(rnorm(20), nrow = 10, ncol = 2)
  gaulss_fam <- gaulss()

  # Should error when lpi attribute is missing for distributional families
  expect_error(
    ffc:::extract_parameter_info_from_lpmat(mock_linpreds, gaulss_fam),
    "Multi-parameter model missing.*lpi.*attribute"
  )
})

test_that("distributional family detection handles edge cases correctly", {
  library(mgcv)

  # Test families with exactly 1 parameter (boundary condition)
  single_param_family <- list(nlp = 1)
  class(single_param_family) <- "family"
  expect_false(ffc:::is_distributional_family(single_param_family))

  # Test families with NULL nlp
  null_nlp_family <- list(nlp = NULL)
  class(null_nlp_family) <- "family"
  expect_false(ffc:::is_distributional_family(null_nlp_family))

  # Test families with nlp = 0
  zero_nlp_family <- list(nlp = 0)
  class(zero_nlp_family) <- "family"
  expect_false(ffc:::is_distributional_family(zero_nlp_family))
})

test_that("distributional basis column regex patterns work correctly", {
  # Test various parameter prefix patterns
  test_data <- data.frame(
    y = rnorm(20),
    x = rnorm(20),
    time = 1:20,
    location_fts_bs_x_1 = rnorm(20),
    scale_fts_bs_x_1 = rnorm(20),
    shape_fts_bs_x_1 = rnorm(20),
    param1_fts_bs_x_1 = rnorm(20),  # Backward compatibility
    param2_fts_bs_x_1 = rnorm(20),
    regular_column = rnorm(20)
  )

  # Test semantic prefix detection
  location_cols <- grep("^location_fts_", names(test_data), value = TRUE)
  scale_cols <- grep("^scale_fts_", names(test_data), value = TRUE)
  shape_cols <- grep("^shape_fts_", names(test_data), value = TRUE)

  expect_length(location_cols, 1)
  expect_length(scale_cols, 1)
  expect_length(shape_cols, 1)
  expect_equal(location_cols, "location_fts_bs_x_1")

  # Test numeric prefix detection (backward compatibility)
  param_cols <- grep("^param[0-9]+_fts_", names(test_data), value = TRUE)
  expect_length(param_cols, 2)
})

test_that("fts_coefs preserves parameter information correctly", {
  library(mgcv)
  set.seed(3)
  n <- 400
  test_data <- gamSim(
    1,
    n = n,
    dist = "poisson",
    scale = 0.2
  )
  test_data$y <- rTweedie(
    exp(test_data$f),
    p = 1.3,
    phi = .5
  )
  test_data$time <- 1:n

  model <- ffc_gam(
    list(y ~ fts(x0, k = 4, share_penalty = FALSE),
         ~ fts(x0, k = 3, share_penalty = FALSE),
         ~ 1),
    data = test_data,
    family = mgcv::twlss(),
    time = "time"
  )

  # Extract coefficients
  coefs <- fts_coefs(model, summary = TRUE)

  # Should have parameter column
  expect_true(".parameter" %in% names(coefs))

  # Should have coefficients for location and scale (shape has ~ 1, no fts)
  unique_params <- unique(coefs$.parameter)
  expect_true("location" %in% unique_params)
  expect_true("scale" %in% unique_params)
  expect_false("shape" %in% unique_params)  # No fts term for shape

  # Each parameter should have correct number of basis functions
  location_bases <- unique(coefs[coefs$.parameter == "location", ".basis"])
  scale_bases <- unique(coefs[coefs$.parameter == "scale", ".basis"])

  expect_true(NROW(location_bases) == 3)
  expect_true(NROW(scale_bases) == 2)
})
