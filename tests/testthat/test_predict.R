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