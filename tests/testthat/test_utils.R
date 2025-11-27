# Tests for utils.R functions - focusing on warnings and errors

test_that("extract_response_vars handles invalid formulas correctly", {
  # Formula with no response
  formula_no_response <- ~ x + y
  
  expect_error(
    ffc:::extract_response_vars(formula_no_response),
    "Formula has no response variable"
  )
  
  # Valid formula should work
  valid_formula <- y ~ x
  result <- ffc:::extract_response_vars(valid_formula)
  expect_equal(result, "y")
})

test_that("extract_base_variables handles empty expressions", {
  # Empty expression should trigger error
  expect_error(
    ffc:::extract_base_variables(""),
    "Empty or whitespace-only expressions not allowed"
  )
  
  # Whitespace only should also trigger error  
  expect_error(
    ffc:::extract_base_variables("   "),
    "Empty or whitespace-only expressions not allowed"
  )
  
  # Valid expression should work
  result <- ffc:::extract_base_variables("log(x)")
  expect_type(result, "character")
})

test_that("compute_functional_predictions handles missing time variable", {
  # Create minimal test data without time variable
  interpreted_data <- data.frame(
    x = 1:5,
    fts_1 = rnorm(5)
  )
  
  functional_fc <- data.frame(
    .time = 1:3,
    .basis = rep(1, 3),
    .sim = rep(1, 3),
    .realisation = rep(1, 3),
    .rep = rep(1, 3),
    fts_1 = rnorm(3)
  )
  
  expect_error(
    ffc:::compute_functional_predictions(
      interpreted_data,
      functional_fc,
      time_var = "missing_time"
    ),
    "Time variable"
  )
})

test_that("compute_functional_predictions handles missing required columns", {
  interpreted_data <- data.frame(
    time = 1:5,
    fts_1 = rnorm(5)
  )
  
  # Missing required columns in functional_fc
  functional_fc_missing <- data.frame(
    .basis = rep(1, 3),
    fts_1 = rnorm(3)
    # Missing .time, .sim, .realisation, .rep
  )
  
  expect_error(
    ffc:::compute_functional_predictions(
      interpreted_data,
      functional_fc_missing,
      time_var = "time"
    ),
    "Missing required columns in functional forecast"
  )
})


test_that("add_row_identifiers handles existing invalid ID column correctly", {
  test_data <- data.frame(
    x = 1:3,
    .row_id = c("a", "a", "c")  # Duplicate values - invalid IDs
  )
  
  # Function should handle invalid IDs by replacing them
  result <- suppressWarnings(ffc:::add_row_identifiers(test_data, ".row_id"))
  
  # Should work and replace the invalid column  
  expect_true(".row_id" %in% names(result))
  expect_equal(nrow(result), 3)
  expect_equal(result$.row_id, 1:3)  # Should be sequential numbers
  
  # Test that valid IDs are preserved
  test_data_valid <- data.frame(
    x = 1:3,
    .row_id = 1:3  # Valid IDs
  )
  
  result_valid <- ffc:::add_row_identifiers(test_data_valid, ".row_id")
  expect_equal(result_valid$.row_id, 1:3)  # Should be unchanged
})

test_that("restore_original_order handles missing ID column", {
  test_data <- data.frame(
    x = 1:3,
    y = 4:6
    # Missing .original_row_id column
  )
  
  expect_error(
    ffc:::restore_original_order(test_data, ".original_row_id"),
    "Row ID column.*original_row_id.*not found"
  )
})

test_that("add_row_identifiers validation mode works", {
  test_data <- data.frame(x = 1:3)
  
  # Validation mode with non-existing column should return FALSE
  expect_silent(
    result <- ffc:::add_row_identifiers(test_data, ".test_id", validate_only = TRUE)
  )
  expect_equal(result, FALSE)  # Column doesn't exist
  
  # Validation mode with valid existing IDs should return TRUE
  test_data_with_valid_id <- data.frame(x = 1:3, .test_id = 1:3)
  expect_silent(
    result <- ffc:::add_row_identifiers(test_data_with_valid_id, ".test_id", validate_only = TRUE)
  )
  expect_equal(result, TRUE)  # Valid IDs exist
})

test_that("utility functions handle edge cases correctly", {
  # Test nlist function
  result <- ffc:::nlist(a = 1, b = 2)
  expect_named(result, c("a", "b"))
  expect_equal(result$a, 1)
  expect_equal(result$b, 2)
  
  # Test c<- function  
  x <- c(1, 2, 3)
  result <- ffc:::`c<-`(x, 4)
  expect_equal(result, c(1, 2, 3, 4))
})

test_that("extract_response_vars handles complex response expressions", {
  # Test cbind response - should return full expression as single character
  formula_cbind <- cbind(y1, y2) ~ x
  result <- ffc:::extract_response_vars(formula_cbind)
  expect_length(result, 1)  # Returns the full cbind expression
  expect_true(grepl("cbind", result))
  
  # Test transformed response
  formula_transform <- log(y) ~ x
  result <- ffc:::extract_response_vars(formula_transform)
  expect_equal(result, "log(y)")
  
  # Test simple response
  formula_simple <- y ~ x
  result <- ffc:::extract_response_vars(formula_simple)
  expect_equal(result, "y")
})