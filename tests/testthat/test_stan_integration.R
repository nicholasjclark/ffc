# Test Stan integration and basic functionality
test_that("Stan model functions are properly exported", {
  # Test that Stan model functions exist
  expect_true(exists("ARDF"))
  expect_true(exists("VARDF"))
  expect_true(exists("GPDF"))
  
  # Test that these are function objects
  expect_true(is.function(ARDF))
  expect_true(is.function(VARDF))
  expect_true(is.function(GPDF))
})

test_that("ARDF function can be called", {
  # Test ARDF function can be invoked without errors
  expect_no_error(ardf_spec <- ARDF(
    formula = .estimate ~ K(K = 2) + p(p = 1),
    family = gaussian(),
    h = 3,
    chains = 1,
    iter = 100
  ))
})

test_that("VARDF function can be called", {
  # Test VARDF function can be invoked without errors
  expect_no_error(vardf_spec <- VARDF(
    formula = .estimate ~ K(K = 2) + p(p = 1),
    family = gaussian(),
    h = 3,
    chains = 1,
    iter = 100
  ))
})

test_that("GPDF function can be called", {
  # Test GPDF function can be invoked without errors
  expect_no_error(gpdf_spec <- GPDF(
    formula = .estimate ~ K(K = 2),
    family = gaussian(),
    h = 3,
    chains = 1,
    iter = 100
  ))
})

test_that("Stan models can accept different parameter values", {
  # Test with different chain numbers
  expect_no_error(ARDF(
    formula = .estimate ~ K(K = 2) + p(p = 1),
    chains = 1
  ))
  
  expect_no_error(ARDF(
    formula = .estimate ~ K(K = 2) + p(p = 1),
    chains = 2
  ))
  
  # Test with different iteration counts
  expect_no_error(ARDF(
    formula = .estimate ~ K(K = 2) + p(p = 1),
    iter = 50
  ))
  
  expect_no_error(ARDF(
    formula = .estimate ~ K(K = 2) + p(p = 1),
    iter = 200
  ))
})

test_that("Stan models handle different horizon specifications", {
  # Test with different forecast horizons
  expect_no_error(ARDF(
    formula = .estimate ~ K(K = 2) + p(p = 1),
    h = 1,
    chains = 1,
    iter = 100
  ))
  
  expect_no_error(ARDF(
    formula = .estimate ~ K(K = 2) + p(p = 1),
    h = 10,
    chains = 1,
    iter = 100
  ))
})

test_that("Stan models can accept different family specifications", {
  # Test with gaussian family
  expect_no_error(ARDF(
    formula = .estimate ~ K(K = 2) + p(p = 1),
    family = gaussian(),
    chains = 1,
    iter = 100
  ))
})

test_that("Stan models handle adaptation parameters", {
  # Test with custom adapt_delta
  expect_no_error(ARDF(
    formula = .estimate ~ K(K = 2) + p(p = 1),
    family = gaussian(),
    adapt_delta = 0.95,
    max_treedepth = 12,
    chains = 1,
    iter = 100
  ))
})

test_that("Stan model training functions exist and have proper structure", {
  # Test that training functions exist
  expect_true(exists("train_ardf", envir = asNamespace("ffc")))
  expect_true(exists("train_vardf", envir = asNamespace("ffc")))
  expect_true(exists("train_gpdf", envir = asNamespace("ffc")))
  
  # Test that these are function objects
  train_ardf_func <- get("train_ardf", envir = asNamespace("ffc"))
  train_vardf_func <- get("train_vardf", envir = asNamespace("ffc"))
  train_gpdf_func <- get("train_gpdf", envir = asNamespace("ffc"))
  
  expect_true(is.function(train_ardf_func))
  expect_true(is.function(train_vardf_func))
  expect_true(is.function(train_gpdf_func))
})

test_that("Stan models have special functions available", {
  # Test ARDF specials exist
  ardf_specials <- ffc:::specials_ardf
  expect_true("K" %in% names(ardf_specials))
  expect_true("p" %in% names(ardf_specials))
  
  # Test VARDF specials exist
  vardf_specials <- ffc:::specials_vardf
  expect_true("K" %in% names(vardf_specials))
  expect_true("p" %in% names(vardf_specials))
  
  # Test GPDF specials exist
  gpdf_specials <- ffc:::specials_gpdf
  expect_true("K" %in% names(gpdf_specials))
})

test_that("forecast integration recognizes Stan model names", {
  # Test that forecast method can handle model name recognition
  functional_coefs <- fts_coefs(example_mod, summary = FALSE, times = 2)
  
  # Test that forecasting function works with standard models
  expect_no_error(SW(forecast(
    functional_coefs,
    model = "ARIMA",
    h = 1,
    times = 2
  )))
  
  # Test model name validation
  expect_error(forecast(
    functional_coefs,
    model = "INVALID_MODEL",
    h = 1,
    times = 2
  ))
})

test_that("Stan model file structure is accessible", {
  # Test that Stan model files can be accessed
  stan_models_dir <- system.file("stan", package = "ffc")
  expect_true(dir.exists(stan_models_dir))
  
  # Test that required stan files exist
  ardf_file <- file.path(stan_models_dir, "ardf.stan")
  vardf_file <- file.path(stan_models_dir, "vardf.stan")  
  gpdf_file <- file.path(stan_models_dir, "gpdf.stan")
  
  expect_true(file.exists(ardf_file))
  expect_true(file.exists(vardf_file))
  expect_true(file.exists(gpdf_file))
})

test_that("prep_tbl_ts_stan works correctly with proper input data", {
  # Create a simple, well-formed test dataset that mimics what prep_tbl_ts_stan expects
  # This avoids the multiple realizations issue by creating clean test data
  set.seed(42)
  n_time <- 10
  basis_names <- c("basis1", "basis2", "basis3")
  
  test_data <- expand.grid(
    .time = 1:n_time,
    .basis = basis_names
  ) %>%
    dplyr::mutate(
      .realisation = 1,  # Single realization to avoid pivot_wider warnings
      .estimate = rnorm(n_time * length(basis_names), 0, 1)
    ) %>%
    tsibble::as_tsibble(
      key = c(.basis, .realisation),
      index = .time
    )
  
  # Test prep_tbl_ts_stan function
  prep_func <- get("prep_tbl_ts_stan", envir = asNamespace("ffc"))
  
  h_value <- 3
  K_value <- 2  # Less than n_series (3)
  p_value <- 1
  
  expect_no_warning({
    stan_data <- prep_func(
      .data = test_data,
      h = h_value,
      K = K_value,
      p = p_value,
      family = gaussian(),
      model = 'ardf'
    )
  })
  
  # Test that required elements are present
  expect_true("n" %in% names(stan_data))
  expect_true("K" %in% names(stan_data))
  expect_true("n_series" %in% names(stan_data))
  expect_true("M" %in% names(stan_data))
  expect_true("n_nonmissing" %in% names(stan_data))
  expect_true("flat_ys" %in% names(stan_data))
  expect_true("obs_ind" %in% names(stan_data))
  expect_true("family" %in% names(stan_data))
  expect_true("alpha" %in% names(stan_data))
  expect_true("P" %in% names(stan_data))
  
  # Test precise dimensions
  n_series_actual <- length(basis_names)
  n_time_expanded <- n_time + h_value
  
  expect_equal(stan_data$n, n_time_expanded)
  expect_equal(stan_data$K, K_value)
  expect_equal(stan_data$n_series, n_series_actual)
  expect_equal(stan_data$P, p_value)
  expect_equal(length(stan_data$alpha), n_series_actual)
  
  # Test M calculation: M = K * (n_series - K) + K * (K - 1) / 2 + K
  expected_M <- K_value * (n_series_actual - K_value) + K_value * (K_value - 1) / 2 + K_value
  expect_equal(stan_data$M, expected_M)
  
  # Test family encoding
  expect_equal(stan_data$family, 1)  # gaussian = 1
  
  # Test data consistency
  expect_equal(length(stan_data$flat_ys), stan_data$n_nonmissing)
  expect_equal(length(stan_data$obs_ind), stan_data$n_nonmissing)
  expect_equal(stan_data$n_nonmissing, n_time * n_series_actual)  # All observations present
  
  # Test that all flat_ys values are numeric and not placeholder values
  expect_true(all(is.numeric(stan_data$flat_ys)))
  expect_true(all(!is.na(stan_data$flat_ys)))
  expect_true(all(stan_data$flat_ys != -1))  # No missing value placeholders
  
  # Test alpha calculation (series means)
  expect_true(all(is.numeric(stan_data$alpha)))
  expect_equal(length(stan_data$alpha), n_series_actual)
  expect_true(all(!is.na(stan_data$alpha)))
})

test_that("prep_tbl_ts_stan handles different model types correctly", {
  # Create simple test data
  set.seed(123)
  test_data <- expand.grid(
    .time = 1:8,
    .basis = c("basis1", "basis2")
  ) %>%
    dplyr::mutate(
      .realisation = 1,
      .estimate = rnorm(16, 0, 1)
    ) %>%
    tsibble::as_tsibble(
      key = c(.basis, .realisation),
      index = .time
    )
  
  prep_func <- get("prep_tbl_ts_stan", envir = asNamespace("ffc"))
  
  # Test ARDF model
  stan_data_ardf <- prep_func(
    test_data, h = 2, K = 1, p = 2,
    family = gaussian(), model = 'ardf'
  )
  expect_true("prior_ar" %in% names(stan_data_ardf))
  expect_true("P" %in% names(stan_data_ardf))
  expect_equal(stan_data_ardf$P, 2)
  
  # Test VARDF model
  stan_data_vardf <- prep_func(
    test_data, h = 2, K = 1, p = 3,
    family = gaussian(), model = 'vardf'
  )
  expect_true("P" %in% names(stan_data_vardf))
  expect_false("prior_ar" %in% names(stan_data_vardf))
  expect_equal(stan_data_vardf$P, 3)
  
  # Test GPDF model
  stan_data_gpdf <- prep_func(
    test_data, h = 2, K = 1, p = 1,
    family = gaussian(), model = 'gpdf'
  )
  expect_false("P" %in% names(stan_data_gpdf))
  expect_false("prior_ar" %in% names(stan_data_gpdf))
})

test_that("prep_tbl_ts_stan validates parameters correctly", {
  # Create minimal test data
  test_data <- expand.grid(
    .time = 1:5,
    .basis = c("basis1")
  ) %>%
    dplyr::mutate(
      .realisation = 1,
      .estimate = rnorm(5)
    ) %>%
    tsibble::as_tsibble(
      key = c(.basis, .realisation),
      index = .time
    )
  
  prep_func <- get("prep_tbl_ts_stan", envir = asNamespace("ffc"))
  
  # Test invalid model name
  expect_error(prep_func(
    test_data, h = 2, K = 1, p = 1,
    family = gaussian(), model = "invalid_model"
  ))
  
  # Test K > n_series warning and adjustment
  expect_warning({
    result <- prep_func(
      test_data, h = 2, K = 5, p = 1,  # K=5 > n_series=1
      family = gaussian(), model = 'ardf'
    )
  }, "K cannot be greater than the number of unique series")
  
  # After warning, K should be adjusted to n_series
  expect_equal(result$K, 1)  # Should be adjusted to n_series
})

test_that("prep_tbl_ts_stan handles missing data correctly", {
  # Create test data with missing values
  set.seed(456)
  test_data <- expand.grid(
    .time = 1:6,
    .basis = c("basis1", "basis2")
  ) %>%
    dplyr::mutate(
      .realisation = 1,
      .estimate = rnorm(12)
    )
  
  # Introduce missing values
  test_data$.estimate[c(2, 8)] <- NA
  
  test_data <- test_data %>%
    tsibble::as_tsibble(
      key = c(.basis, .realisation),
      index = .time
    )
  
  prep_func <- get("prep_tbl_ts_stan", envir = asNamespace("ffc"))
  
  stan_data <- prep_func(
    test_data, h = 2, K = 1, p = 1,
    family = gaussian(), model = 'ardf'
  )
  
  # Test that missing data is handled correctly
  total_possible_obs <- 6 * 2  # n_time * n_series
  expect_true(stan_data$n_nonmissing < total_possible_obs)
  expect_equal(stan_data$n_nonmissing, total_possible_obs - 2)  # 2 missing values
  
  # Test that flat_ys doesn't contain missing value placeholders
  expect_equal(length(stan_data$flat_ys), stan_data$n_nonmissing)
  expect_true(all(!is.na(stan_data$flat_ys)))
  expect_true(all(stan_data$flat_ys != -1))
})

test_that("Stan data preparation functions exist and can be called", {
  # Test that internal Stan data preparation functions exist
  expect_true(exists("prep_tbl_ts_stan", envir = asNamespace("ffc")))
  expect_true(exists("extract_stan_fc", envir = asNamespace("ffc")))
  
  # Test that these are function objects
  prep_func <- get("prep_tbl_ts_stan", envir = asNamespace("ffc"))
  extract_func <- get("extract_stan_fc", envir = asNamespace("ffc"))
  
  expect_true(is.function(prep_func))
  expect_true(is.function(extract_func))
  
  # Test that prep function has expected parameters
  prep_formals <- names(formals(prep_func))
  expected_params <- c(".data", "h", "K", "p", "family", "model")
  expect_true(all(expected_params %in% prep_formals))
})

test_that("binomial models work with cbind() response", {
  # Create simulated binomial data with trials and successes
  set.seed(101112)
  n_obs <- 50
  
  binomial_data <- data.frame(
    time = 1:n_obs,
    x = rnorm(n_obs),
    season = rep(1:4, length.out = n_obs),
    trials = sample(20:100, n_obs, replace = TRUE)
  )
  
  # Generate successes based on logistic model
  linear_pred <- -1 + 0.5 * binomial_data$x + 
    0.3 * sin(2 * pi * binomial_data$season / 4)
  prob <- plogis(linear_pred)
  binomial_data$successes <- rbinom(n_obs, binomial_data$trials, prob)
  
  # Add failures column for proper cbind syntax
  binomial_data$failures <- binomial_data$trials - binomial_data$successes
  
  # Test model fitting with cbind() response
  expect_no_error({
    mod_cbind <- ffc_gam(
      cbind(successes, failures) ~ 
        fts(time, mean_only = TRUE, time_k = 10) +
        fts(season, bs = 'cc', k = 4),
      knots = list(season = c(0.5, 4.5)),
      data = binomial_data,
      time = "time",
      family = binomial()
    )
  })
  
  # Test that model structure is correct
  expect_true(inherits(mod_cbind, "ffc_gam"))
  expect_true(inherits(mod_cbind, "gam"))
  expect_true(mod_cbind$family$family == "binomial")
  
  # Test predictions work
  expect_no_error({
    preds_link <- predict(mod_cbind, type = "link")
    preds_response <- predict(mod_cbind, type = "response")
  })
  
  expect_true(is.numeric(preds_link))
  expect_true(is.numeric(preds_response))
  expect_true(all(preds_response >= 0 & preds_response <= 1))
  
  # Test model frame extraction
  expect_no_error({
    mod_frame <- model.frame(mod_cbind)
  })
  
  expect_true(is.data.frame(mod_frame))
  
  # Test coefficient extraction
  expect_no_error({
    fts_coefs_result <- fts_coefs(mod_cbind, summary = TRUE)
  })
  
  expect_true(inherits(fts_coefs_result, "fts_ts"))
})

test_that("binomial models handle response parsing correctly", {
  # Test that cbind response parsing works correctly
  set.seed(131415)
  n_obs <- 30
  
  test_data <- data.frame(
    time = 1:n_obs,
    var1 = rnorm(n_obs),
    trials = sample(10:50, n_obs, replace = TRUE)
  )
  test_data$successes <- rbinom(n_obs, test_data$trials, 0.3)
  test_data$failures <- test_data$trials - test_data$successes
  
  # Fit model with cbind response
  mod <- ffc_gam(
    cbind(successes, failures) ~ 
      fts(time, mean_only = TRUE, time_k = 8),
    data = test_data,
    time = "time",
    family = binomial()
  )
  
  # Test that response variables are correctly identified
  # The function should handle cbind parsing internally
  expect_true(inherits(mod, "ffc_gam"))
  
  # Test that missing data handling works with cbind
  test_data_na <- test_data
  test_data_na$successes[c(5, 15)] <- NA
  test_data_na$failures <- test_data_na$trials - test_data_na$successes
  
  expect_no_error({
    mod_na <- ffc_gam(
      cbind(successes, failures) ~ 
        fts(time, mean_only = TRUE, time_k = 8),
      data = test_data_na,
      time = "time",
      family = binomial()
    )
  })
  
  expect_true(inherits(mod_na, "ffc_gam"))
})