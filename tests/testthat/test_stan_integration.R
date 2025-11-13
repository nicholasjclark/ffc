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