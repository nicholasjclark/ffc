# Test prediction functionality
test_that("predict.ffc_gam() works without newdata", {
  # Test basic prediction without newdata (uses training data)
  pred_link <- predict(example_mod, type = "link")
  pred_response <- predict(example_mod, type = "response")
  pred_terms <- predict(example_mod, type = "terms")
  
  # Check basic structure
  expect_true(is.numeric(pred_link))
  expect_true(is.numeric(pred_response))
  expect_true(is.matrix(pred_terms))
  
  # Check dimensions match training data
  expect_equal(length(pred_link), nrow(example_mod$model))
  expect_equal(length(pred_response), nrow(example_mod$model))
  expect_equal(nrow(pred_terms), nrow(example_mod$model))
  
  # Response predictions should be non-negative for Poisson
  expect_range(pred_response, lower = 0)
})

test_that("predict.ffc_gam() works with newdata", {
  # Create test newdata
  test_newdata <- data.frame(
    season = c(1, 6, 12),
    time = c(76, 77, 78)
  )
  
  # Test predictions with newdata
  pred_link <- predict(example_mod, newdata = test_newdata, type = "link")
  pred_response <- predict(example_mod, newdata = test_newdata, type = "response")
  
  # Check structure
  expect_true(is.numeric(pred_link))
  expect_true(is.numeric(pred_response))
  expect_equal(length(pred_link), nrow(test_newdata))
  expect_equal(length(pred_response), nrow(test_newdata))
  
  # Response predictions should be non-negative for Poisson
  expect_range(pred_response, lower = 0)
})

test_that("predict.ffc_gam() handles different prediction types", {
  test_newdata <- data.frame(
    season = c(1, 6),
    time = c(76, 77)
  )
  
  # Test all prediction types
  pred_link <- predict(example_mod, newdata = test_newdata, type = "link")
  pred_response <- predict(example_mod, newdata = test_newdata, type = "response") 
  pred_terms <- predict(example_mod, newdata = test_newdata, type = "terms")
  pred_lpmatrix <- predict(example_mod, newdata = test_newdata, type = "lpmatrix")
  pred_iterms <- predict(example_mod, newdata = test_newdata, type = "iterms")
  
  # Check types and dimensions
  expect_true(is.numeric(pred_link))
  expect_true(is.numeric(pred_response))
  expect_true(is.matrix(pred_terms))
  expect_true(is.matrix(pred_lpmatrix))
  expect_true(is.matrix(pred_iterms))
  
  # Check dimensions
  expect_equal(length(pred_link), nrow(test_newdata))
  expect_equal(length(pred_response), nrow(test_newdata))
  expect_equal(nrow(pred_terms), nrow(test_newdata))
  expect_equal(nrow(pred_lpmatrix), nrow(test_newdata))
  expect_equal(nrow(pred_iterms), nrow(test_newdata))
  
  # lpmatrix should have columns for each coefficient
  expect_equal(ncol(pred_lpmatrix), length(coef(example_mod)))
})

test_that("predict.ffc_gam() works with se.fit = TRUE", {
  test_newdata <- data.frame(
    season = c(1, 6),
    time = c(76, 77)
  )
  
  # Test with standard errors
  pred_se_link <- predict(example_mod, newdata = test_newdata, 
                          type = "link", se.fit = TRUE)
  pred_se_response <- predict(example_mod, newdata = test_newdata,
                              type = "response", se.fit = TRUE)
  
  # Should return list with fit and se.fit components
  expect_true(is.list(pred_se_link))
  expect_true(is.list(pred_se_response))
  expect_true(all(c("fit", "se.fit") %in% names(pred_se_link)))
  expect_true(all(c("fit", "se.fit") %in% names(pred_se_response)))
  
  # Check dimensions
  expect_equal(length(pred_se_link$fit), nrow(test_newdata))
  expect_equal(length(pred_se_link$se.fit), nrow(test_newdata))
  expect_equal(length(pred_se_response$fit), nrow(test_newdata))
  expect_equal(length(pred_se_response$se.fit), nrow(test_newdata))
  
  # Standard errors should be non-negative
  expect_range(pred_se_link$se.fit, lower = 0)
  expect_range(pred_se_response$se.fit, lower = 0)
})

test_that("predict.ffc_gam() validates type parameter", {
  test_newdata <- data.frame(
    season = c(1),
    time = c(76)
  )
  
  # Test invalid type
  expect_error(predict(example_mod, newdata = test_newdata, type = "invalid"))
})

test_that("predict.ffc_gam() handles missing variables in newdata", {
  # Test missing season variable
  incomplete_newdata1 <- data.frame(time = c(76, 77))
  expect_error(SW(predict(example_mod, newdata = incomplete_newdata1)))
  
  # Test missing time variable  
  incomplete_newdata2 <- data.frame(season = c(1, 6))
  expect_error(SW(predict(example_mod, newdata = incomplete_newdata2)))
})

test_that("predict.ffc_gam() handles edge cases in newdata", {
  # Test single observation
  single_obs <- data.frame(season = 1, time = 76)
  pred_single <- predict(example_mod, newdata = single_obs, type = "response")
  expect_equal(length(pred_single), 1)
  expect_true(is.numeric(pred_single))
  expect_range(pred_single, lower = 0)
  
  # Test with season values that match training data range
  extended_newdata <- data.frame(
    season = c(2, 7, 11),  # Use numeric values within training range
    time = c(76, 77, 78)
  )
  pred_extended <- predict(example_mod, newdata = extended_newdata, type = "response")
  expect_equal(length(pred_extended), nrow(extended_newdata))
  expect_range(pred_extended, lower = 0)
})

test_that("predict.ffc_gam() handles newdata with different time ranges", {
  # Test predictions beyond training time range
  future_newdata <- data.frame(
    season = c(1, 6, 12),
    time = c(100, 101, 102)  # Well beyond training range
  )
  pred_future <- predict(example_mod, newdata = future_newdata, type = "response")
  expect_equal(length(pred_future), nrow(future_newdata))
  expect_true(all(is.finite(pred_future)))
  expect_range(pred_future, lower = 0)
  
  # Test predictions before training time range  
  past_newdata <- data.frame(
    season = c(1, 6, 12),
    time = c(-5, -4, -3)  # Before training range
  )
  pred_past <- predict(example_mod, newdata = past_newdata, type = "response")
  expect_equal(length(pred_past), nrow(past_newdata))
  expect_true(all(is.finite(pred_past)))
  expect_range(pred_past, lower = 0)
})

test_that("predict.ffc_gam() preserves prediction consistency", {
  test_newdata <- data.frame(
    season = c(1, 6, 12),
    time = c(76, 77, 78)
  )
  
  # Same newdata should give same predictions
  pred1 <- predict(example_mod, newdata = test_newdata, type = "link")
  pred2 <- predict(example_mod, newdata = test_newdata, type = "link")
  expect_equal(pred1, pred2)
  
  # Different order should still give same predictions (when reordered)
  reordered_newdata <- test_newdata[c(3, 1, 2), ]
  pred_reordered <- predict(example_mod, newdata = reordered_newdata, type = "link")
  expect_equal(pred_reordered[c(2, 3, 1)], pred1)
})

test_that("predict.ffc_gam() works with terms predictions", {
  test_newdata <- data.frame(
    season = c(1, 6, 12),
    time = c(76, 77, 78)
  )
  
  pred_terms <- predict(example_mod, newdata = test_newdata, type = "terms")
  
  # Should be a matrix
  expect_true(is.matrix(pred_terms))
  expect_equal(nrow(pred_terms), nrow(test_newdata))
  
  # Should have columns for each smooth term (but not intercept)
  expect_true(ncol(pred_terms) > 0)
  
  # Column names should indicate the terms
  expect_true(!is.null(colnames(pred_terms)))
})

test_that("predict.ffc_gam() handles lpmatrix correctly", {
  test_newdata <- data.frame(
    season = c(1, 6),
    time = c(76, 77)  
  )
  
  pred_lpmatrix <- predict(example_mod, newdata = test_newdata, type = "lpmatrix")
  
  # Should be a matrix with rows = newdata, cols = coefficients
  expect_true(is.matrix(pred_lpmatrix))
  expect_equal(nrow(pred_lpmatrix), nrow(test_newdata))
  expect_equal(ncol(pred_lpmatrix), length(coef(example_mod)))
  
  # Should have model.offset attribute if offset was used
  expect_true(!is.null(attr(pred_lpmatrix, "model.offset")))
  
  # Matrix multiplication with coefficients should give link predictions
  manual_pred <- as.vector(pred_lpmatrix %*% coef(example_mod) + 
                           attr(pred_lpmatrix, "model.offset"))
  auto_pred <- as.vector(predict(example_mod, newdata = test_newdata, type = "link"))
  expect_equal(manual_pred, auto_pred, tolerance = 1e-10)
})

test_that("predict.ffc_gam() integrates correctly with fts basis functions", {
  # Create newdata that will trigger fts basis function evaluation
  test_newdata <- data.frame(
    season = c(1, 6, 12),
    time = c(76, 77, 78)  # Future time points
  )
  
  # Predictions should work without error
  pred <- predict(example_mod, newdata = test_newdata, type = "response")
  expect_equal(length(pred), nrow(test_newdata))
  expect_true(all(is.finite(pred)))
  expect_range(pred, lower = 0)
  
  # lpmatrix should include fts basis evaluations
  lpmat <- predict(example_mod, newdata = test_newdata, type = "lpmatrix") 
  expect_true(is.matrix(lpmat))
  expect_equal(nrow(lpmat), nrow(test_newdata))
})

# Additional distributional prediction tests
test_that("split_linear_predictors_by_lpi validates parameter indices", {
  library(mgcv)
  
  # Create linpreds matrix with 10 columns
  linpreds <- matrix(rnorm(100), nrow = 10, ncol = 10)
  
  # Create invalid parameter info with out-of-bounds indices
  invalid_param_info <- list(
    n_parameters = 2L,
    parameter_indices = list(1:5, 11:15)  # Second set exceeds ncol(linpreds)
  )
  
  # Should error when indices exceed matrix columns
  expect_error(
    ffc:::split_linear_predictors_by_lpi(linpreds, invalid_param_info),
    "Assertion on 'indices' failed"
  )
})

test_that("apply_distributional_inverse_links validates parameter count", {
  library(mgcv)
  
  # Create gaulss family (2 parameters)
  gaulss_fam <- gaulss()
  
  # Create parameter predictions with wrong number of parameters
  wrong_par_predictions <- list(
    location = matrix(rnorm(10), nrow = 5),
    scale = matrix(rnorm(10), nrow = 5), 
    shape = matrix(rnorm(10), nrow = 5)   # gaulss only has 2 parameters
  )
  
  # Should error when parameter count exceeds family nlp
  expect_error(
    ffc:::apply_distributional_inverse_links(wrong_par_predictions, gaulss_fam),
    "par_predictions length vs family parameters"
  )
})

test_that("distributional families handle invalid link function calls", {
  library(mgcv)
  
  # Create mock parameter predictions with invalid dimensions
  invalid_predictions <- list(
    location = matrix(c(Inf, -Inf, NaN), nrow = 3),  # Invalid values
    scale = matrix(c(1, 2, 3), nrow = 3)
  )
  
  gaulss_fam <- gaulss()
  
  # Should handle infinite/NaN values gracefully
  result <- SW(ffc:::apply_distributional_inverse_links(invalid_predictions, gaulss_fam))
  
  # Result should be a list with same structure
  expect_true(is.list(result))
  expect_length(result, 2)
  
  # Test with completely empty predictions
  empty_predictions <- list(
    location = matrix(numeric(0), nrow = 0),
    scale = matrix(numeric(0), nrow = 0)
  )
  
  empty_result <- ffc:::apply_distributional_inverse_links(empty_predictions, gaulss_fam)
  expect_true(is.list(empty_result))
  expect_length(empty_result, 2)
  expect_equal(nrow(empty_result[[1]]), 0)
  expect_equal(nrow(empty_result[[2]]), 0)
})

# Tests for new posterior_predict parameter matrix functionality
test_that("posterior_predict validates parameter matrix dimensions", {
  set.seed(42)
  n_obs <- 5
  linpreds <- matrix(rnorm(n_obs * 2), nrow = n_obs, ncol = 2)
  
  # Test wrong row dimensions should error before function execution
  wrong_matrix <- matrix(rnorm(3 * 2), nrow = 3, ncol = 2)
  
  # This should error in validation step (checkmate assertion)
  expect_error({
    ffc:::posterior_predict(
      object = list(family = gaussian()),
      linpreds = linpreds,
      location_matrix = wrong_matrix
    )
  }, "Must have exactly 5 rows")
})

test_that("posterior_predict validates parameter matrix column consistency", {
  set.seed(42)
  n_obs <- 5
  linpreds <- matrix(rnorm(n_obs * 2), nrow = n_obs, ncol = 2)
  
  # Test inconsistent column counts between parameter matrices
  location_matrix <- matrix(rnorm(n_obs * 3), nrow = n_obs, ncol = 3)
  scale_matrix <- matrix(rnorm(n_obs * 2), nrow = n_obs, ncol = 2)
  
  expect_error({
    ffc:::posterior_predict(
      object = list(family = gaussian()),
      linpreds = linpreds,
      location_matrix = location_matrix,
      scale_matrix = scale_matrix
    )
  }, "same number of columns")
})

test_that("posterior_predict validates non-matrix parameter inputs", {
  set.seed(42)
  n_obs <- 5
  linpreds <- matrix(rnorm(n_obs * 2), nrow = n_obs, ncol = 2)
  
  # Non-matrix parameter should error in validation
  expect_error({
    ffc:::posterior_predict(
      object = list(family = gaussian()),
      linpreds = linpreds,
      location_matrix = c(1, 2, 3, 4, 5)  # vector not matrix
    )
  }, "matrix")
})

test_that("posterior_predict validates family-specific parameter requirements", {
  library(mgcv)
  set.seed(42)
  n_obs <- 5
  linpreds <- matrix(rnorm(n_obs * 2), nrow = n_obs, ncol = 2)
  
  # Test gaulss family (nlp=2, requires location + scale)
  gaulss_family <- gaulss()
  gaulss_object <- list(family = gaulss_family)
  
  location_matrix <- matrix(rnorm(n_obs * 3), nrow = n_obs, ncol = 3)
  scale_matrix <- matrix(rnorm(n_obs * 3), nrow = n_obs, ncol = 3)
  shape_matrix <- matrix(rnorm(n_obs * 3), nrow = n_obs, ncol = 3)
  
  # Missing scale_matrix should error
  expect_error({
    ffc:::posterior_predict(gaulss_object, linpreds,
                           location_matrix = location_matrix)
  }, "gaulss.*requires.*scale_matrix")
  
  # Extra shape_matrix should error  
  expect_error({
    ffc:::posterior_predict(gaulss_object, linpreds,
                           location_matrix = location_matrix,
                           scale_matrix = scale_matrix,
                           shape_matrix = shape_matrix)
  }, "gaulss.*does not use.*shape_matrix")
})

test_that("posterior_predict validates twlss family parameter requirements", {
  library(mgcv)
  set.seed(42)
  n_obs <- 5
  linpreds <- matrix(rnorm(n_obs * 3), nrow = n_obs, ncol = 3)
  
  # Test twlss family (nlp=3, requires location + scale + shape)
  twlss_family <- twlss()
  twlss_object <- list(family = twlss_family)
  
  location_matrix <- matrix(rnorm(n_obs * 2), nrow = n_obs, ncol = 2)
  scale_matrix <- matrix(rnorm(n_obs * 2), nrow = n_obs, ncol = 2)
  
  # Missing shape_matrix should error
  expect_error({
    ffc:::posterior_predict(twlss_object, linpreds,
                           location_matrix = location_matrix,
                           scale_matrix = scale_matrix)
  }, "twlss.*requires.*shape_matrix")
})

