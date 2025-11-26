# Test data handling and edge cases
test_that("interpret_ffc() handles tibble to data.frame conversion", {
  # Test with regular data.frame (should work)
  dat_df <- as.data.frame(dat)
  expect_no_error(interpret_ffc(
    formula = y ~ s(season, bs = "cc", k = 12) + fts(time, time_k = 5),
    data = dat_df,
    time_var = "time"
  ))

  # Test with tibble (should convert to data.frame internally)
  dat_tbl <- tibble::as_tibble(dat)
  result_tbl <- interpret_ffc(
    formula = y ~ s(season, bs = "cc", k = 12) + fts(time, time_k = 5),
    data = dat_tbl,
    time_var = "time"
  )

  # Should work without error
  expect_true(inherits(result_tbl, "list"))
  expect_true("data" %in% names(result_tbl))
  expect_true(inherits(result_tbl$data, "data.frame"))
})

test_that("interpret_ffc() preserves original data structure", {
  dat_tbl <- tibble::as_tibble(dat)
  result <- interpret_ffc(
    formula = y ~ s(season, bs = "cc", k = 12) + fts(time, time_k = 5),
    data = dat_tbl,
    time_var = "time"
  )

  # Original data should be preserved
  expect_true("orig_data" %in% names(result))
  expect_true(inherits(result$orig_data, "tbl_df"))
  expect_equal(nrow(result$orig_data), nrow(dat_tbl))
  expect_equal(ncol(result$orig_data), ncol(dat_tbl))
})

test_that("expand_tbl_ts() works correctly", {
  # Create a test fts_ts object
  test_coefs <- fts_coefs(example_mod, summary = FALSE, times = 2)

  # Convert to the format expected by expand_tbl_ts
  test_data <- test_coefs %>%
    dplyr::select(.basis, .realisation, .time = time, .estimate)

  # Test expansion
  expanded <- expand_tbl_ts(test_data, h = 3)

  # Should have more rows due to expansion
  expect_true(nrow(expanded) > nrow(test_data))

  # Should include NA values for future time points
  expect_true(any(is.na(expanded$.estimate)))

  # Time range should be extended
  original_max_time <- max(test_data$.time)
  expanded_max_time <- max(expanded$.time)
  expect_equal(expanded_max_time, original_max_time + 3)
})

test_that("make_future_data() creates appropriate future data structure", {
  # Create a simple tsibble for testing
  test_tsibble <- tsibble::tsibble(
    time = 1:10,
    value = rnorm(10),
    index = time
  )

  # Test future data creation
  future_data <- make_future_data(test_tsibble, h = 5)

  # Should be a tsibble
  expect_true(tsibble::is_tsibble(future_data))

  # Should have 5 future time points
  expect_equal(nrow(future_data), 5)

  # Time should continue from where original data ended
  expect_true(min(future_data$time) > max(test_tsibble$time))
})

test_that("tsibble checking functions work correctly", {
  # Test check_gaps with data that has no gaps
  regular_data <- tsibble::tsibble(
    time = 1:10,
    value = rnorm(10),
    index = time
  )
  expect_no_error(check_gaps(regular_data))

  # Test check_regular with regular data
  expect_no_error(check_regular(regular_data))

  # Test check_ordered with ordered data
  expect_no_error(check_ordered(regular_data))

  # Test all_tsbl_checks with valid data
  expect_no_error(all_tsbl_checks(regular_data))
})

test_that("tsibble checking functions detect issues", {
  # Test with irregular data
  irregular_data <- tsibble::tsibble(
    time = c(1, 3, 7, 10),  # irregular spacing
    value = rnorm(4),
    index = time,
    regular = FALSE
  )
  expect_error(check_regular(irregular_data), "irregular time series")

  # Test with empty data
  empty_data <- tsibble::tsibble(
    time = integer(0),
    value = numeric(0),
    index = time
  )
  expect_error(all_tsbl_checks(empty_data), "no data to model")
})

test_that("ffc_gam() handles different data types", {
  # Test with data.frame
  dat_df <- as.data.frame(dat)
  mod_df <- ffc_gam(
    y ~ s(season, bs = "cc", k = 12) + fts(time, time_k = 5),
    time = "time",
    data = dat_df,
    family = poisson()
  )
  expect_true(inherits(mod_df, "ffc_gam"))

  # Test with tibble
  dat_tbl <- tibble::as_tibble(dat)
  mod_tbl <- ffc_gam(
    y ~ s(season, bs = "cc", k = 12) + fts(time, time_k = 5),
    time = "time",
    data = dat_tbl,
    family = poisson()
  )
  expect_true(inherits(mod_tbl, "ffc_gam"))

  # Results should be equivalent
  expect_equal(length(coef(mod_df)), length(coef(mod_tbl)))
  expect_equal(nrow(mod_df$model), nrow(mod_tbl$model))
})

test_that("ffc_gam() validates time variable exists", {
  # Test with missing time variable
  expect_error(
    ffc_gam(
      y ~ s(season, bs = "cc", k = 12) + fts(time, time_k = 5),
      time = "nonexistent_time",
      data = dat,
      family = poisson()
    )
  )
})

test_that("missing data handling works correctly", {
  # Create data with missing values
  dat_missing <- dat
  dat_missing$y[c(10, 20, 30)] <- NA
  dat_missing$season[c(15, 25)] <- NA

  # Should reject missing values with helpful error
  expect_error(
    ffc_gam(
      y ~ s(season, bs = "cc", k = 12) + fts(time, time_k = 5),
      time = "time",
      data = dat_missing,
      family = poisson()
    ),
    "Missing values detected in variables"
  )
})

test_that("offset handling works correctly", {
  # Create data with offset variable
  dat_offset <- dat
  dat_offset$offset_var <- log(dat_offset$y + 1)

  # Test with offset in formula
  mod_offset <- ffc_gam(
    y ~ offset(offset_var) + s(season, bs = "cc", k = 12) + fts(time, time_k = 5),
    time = "time",
    data = dat_offset,
    family = poisson()
  )

  expect_true(inherits(mod_offset, "ffc_gam"))
  expect_true(!is.null(mod_offset$offset))
})

test_that("formula interpretation handles complex cases", {
  # Test formula with single fts term (avoid multiple smooths causing df issues)
  complex_result <- interpret_ffc(
    formula = y ~ s(season, k = 5) + fts(time, time_k = 5),
    data = dat,
    time_var = "time"
  )

  expect_true(inherits(complex_result, "list"))
  expect_true(length(complex_result$fts_smooths) == 1)

  # Test formula with no fts terms
  simple_result <- interpret_ffc(
    formula = y ~ s(season, k = 5) + s(year, k = 4),
    data = dat,
    time_var = "time"
  )

  expect_true(inherits(simple_result, "list"))
  expect_true(length(simple_result$fts_smooths) == 0)
})

test_that("edge cases in data dimensions", {
  # Test with larger minimal data to avoid degrees of freedom issues
  minimal_dat <- data.frame(
    y = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
    season = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
    time = 1:10
  )

  # Should work but may produce warnings about small data size
  expect_no_error(
    SW(mod_minimal <- ffc_gam(
      y ~ s(season, k = 3) + fts(time, time_k = 3),
      time = "time",
      data = minimal_dat,
      family = poisson()
    ))
  )

  expect_true(inherits(mod_minimal, "ffc_gam"))
})

test_that("data structure preservation in predictions", {
  # Test that newdata structure is handled correctly in predictions
  newdata_df <- data.frame(
    season = c(1, 6, 12),
    time = c(76, 77, 78)
  )

  newdata_tbl <- tibble::as_tibble(newdata_df)

  # Predictions should work with both data types
  pred_df <- predict(example_mod, newdata = newdata_df, type = "response")
  pred_tbl <- predict(example_mod, newdata = newdata_tbl, type = "response")

  # Results should be identical
  expect_equal(pred_df, pred_tbl)
  expect_equal(length(pred_df), nrow(newdata_df))
  expect_equal(length(pred_tbl), nrow(newdata_tbl))
})

test_that("fts_coefs() handles different summary options", {
  # Test summarized coefficients
  coefs_summary <- fts_coefs(example_mod, summary = TRUE)
  expect_true(inherits(coefs_summary, "fts_ts"))
  expect_true(attr(coefs_summary, "summarized"))
  expect_true(all(c(".estimate", ".se") %in% names(coefs_summary)))

  # Test raw coefficients
  coefs_raw <- fts_coefs(example_mod, summary = FALSE, n_samples = 3)
  expect_true(inherits(coefs_raw, "fts_ts"))
  expect_false(attr(coefs_raw, "summarized"))
  expect_true(".realisation" %in% names(coefs_raw))
  expect_true(length(unique(coefs_raw$.realisation)) == 3)
})

test_that("model.frame.ffc_gam() works correctly", {
  # Test model frame extraction
  mf <- model.frame(example_mod)

  expect_true(inherits(mf, "data.frame"))
  expect_equal(nrow(mf), nrow(example_mod$model))
  expect_true(all(names(example_mod$model) %in% names(mf)))
})

test_that("Stan data preparation function exists and has proper structure", {
  # Test that prep_tbl_ts_stan function exists and can be called
  expect_true(exists("prep_tbl_ts_stan", envir = asNamespace("ffc")))

  # Test that function has proper parameters
  prep_func <- get("prep_tbl_ts_stan", envir = asNamespace("ffc"))
  expect_true(is.function(prep_func))

  # Test argument validation
  expect_error(prep_tbl_ts_stan(
    .data = data.frame(),
    h = 2,
    K = 1,
    p = 1,
    family = gaussian(),
    model = "invalid_model"
  ))
})
