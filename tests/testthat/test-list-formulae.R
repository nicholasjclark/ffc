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