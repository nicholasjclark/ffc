# Test input validation integration in ffc package functions
test_that("fts() validates time_k parameter correctly", {
  # Test valid time_k parameter
  expect_no_error(fts(time, time_k = 1))
  expect_no_error(fts(time, time_k = 5))
  expect_no_error(fts(time, time_k = 20))
  
  # Test invalid time_k parameter
  expect_error(fts(time, time_k = 0))
  expect_error(fts(time, time_k = -1)) 
  expect_error(fts(time, time_k = 1.5))
  expect_error(fts(time, time_k = "5"))
  expect_error(fts(time, time_k = NA))
  expect_error(fts(time, time_k = NULL))
})

test_that("forecast.fts_ts() validates h and times parameters correctly", {
  functional_coefs <- fts_coefs(example_mod, summary = FALSE, times = 2)
  
  # Test valid h and times parameters
  expect_no_error(SW(forecast(functional_coefs, h = 1, times = 1, model = "ARIMA")))
  expect_no_error(SW(forecast(functional_coefs, h = 5, times = 10, model = "ARIMA")))
  
  # Test invalid h parameter
  expect_error(forecast(functional_coefs, h = 0, times = 2, model = "ARIMA"))
  expect_error(forecast(functional_coefs, h = -1, times = 2, model = "ARIMA"))
  expect_error(forecast(functional_coefs, h = 1.5, times = 2, model = "ARIMA"))
  expect_error(forecast(functional_coefs, h = "2", times = 2, model = "ARIMA"))
  
  # Test invalid times parameter
  expect_error(forecast(functional_coefs, h = 2, times = 0, model = "ARIMA"))
  expect_error(forecast(functional_coefs, h = 2, times = -1, model = "ARIMA"))
  expect_error(forecast(functional_coefs, h = 2, times = 1.5, model = "ARIMA"))
  expect_error(forecast(functional_coefs, h = 2, times = "5", model = "ARIMA"))
})

test_that("fts_coefs() validates times parameter correctly", {
  # Test valid times parameter (use >= 2 to avoid matrix dimension issues)
  expect_no_error(fts_coefs(example_mod, summary = FALSE, times = 2))
  expect_no_error(fts_coefs(example_mod, summary = FALSE, times = 5))
  expect_no_error(fts_coefs(example_mod, summary = FALSE, times = 100))
  
  # Test invalid times parameter
  expect_error(fts_coefs(example_mod, summary = FALSE, times = 0))
  expect_error(fts_coefs(example_mod, summary = FALSE, times = -1))
  expect_error(fts_coefs(example_mod, summary = FALSE, times = 1.5))
  expect_error(fts_coefs(example_mod, summary = FALSE, times = "5"))
  expect_error(fts_coefs(example_mod, summary = FALSE, times = NA))
})

test_that("validation provides clear error messages", {
  # Test that validation errors contain helpful information
  expect_error(fts(time, time_k = 0), "time_k.*>= 1")
  expect_error(fts(time, time_k = -1), "time_k.*>= 1")
  
  functional_coefs <- fts_coefs(example_mod, summary = FALSE, times = 2)
  expect_error(forecast(functional_coefs, h = 0, model = "ARIMA"), "h.*>= 1")
  expect_error(forecast(functional_coefs, times = 0, model = "ARIMA"), "times.*>= 1")
  
  expect_error(fts_coefs(example_mod, summary = FALSE, times = 0), "times.*>= 1")
})

test_that("validation handles boundary cases correctly", {
  # Test minimum valid values
  expect_no_error(fts(time, time_k = 1))
  
  functional_coefs <- fts_coefs(example_mod, summary = FALSE, times = 2)
  expect_no_error(SW(forecast(functional_coefs, h = 1, times = 1, model = "ARIMA")))
  
  expect_no_error(fts_coefs(example_mod, summary = FALSE, times = 2))
  
  # Test reasonably large values
  expect_no_error(fts(time, time_k = 50))
  expect_no_error(SW(forecast(functional_coefs, h = 10, times = 50, model = "ARIMA")))
  expect_no_error(fts_coefs(example_mod, summary = FALSE, times = 50))
})

test_that("validation integrates properly with function flow", {
  # Test that validation doesn't interfere with normal function operation
  
  # fts() with valid parameters should work normally
  fts_call <- fts(time, time_k = 5)
  expect_true(inherits(fts_call, "list"))
  
  # forecast() with valid parameters should work normally
  functional_coefs <- fts_coefs(example_mod, summary = FALSE, times = 3)
  fc <- SW(forecast(functional_coefs, h = 2, times = 3, model = "ARIMA"))
  expect_true(inherits(fc, "tbl_df"))
  
  # fts_coefs() with valid parameters should work normally
  coefs <- fts_coefs(example_mod, summary = FALSE, times = 3)
  expect_true(inherits(coefs, "fts_ts"))
})