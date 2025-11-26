# Tests for Stan-based dynamic factor models (ARDF, GPDF, VARDF)
# These tests are skipped on CRAN due to computational intensity

test_that("ARDF model works with basic functionality", {
  skip_on_cran()

  # Simple test data
  test_data <- data.frame(
    y = c(
      10,
      12,
      14,
      13,
      15,
      16,
      18,
      17,
      19,
      20,
      22,
      21,
      23,
      24,
      26,
      25,
      27,
      28,
      30,
      29
    ),
    season = rep(1:4, 5),
    time = 1:20
  )

  # Fit simple model
  mod <- SW(ffc_gam(
    y ~ fts(time, mean_only = TRUE, time_k = 8, k = 3) +
      fts(season, k = 3, time_k = 8),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  # Test ARDF forecasting
  newdata <- data.frame(
    season = c(1, 2),
    time = c(21, 22)
  )

  # Basic ARDF forecast
  fc_ardf <- SW(forecast(
    mod,
    newdata = newdata,
    model = "ARDF",
    K = 1,
    lag = 1,
    chains = 1,
    iter = 300,
    summary = TRUE
  ))

  fc_ardf$.estimate
  fc_ardf$.error

  expect_true(inherits(fc_ardf, "tbl_df"))
  expect_true(all(c(".estimate", ".error") %in% names(fc_ardf)))
  expect_equal(nrow(fc_ardf), nrow(newdata))
  expect_true(all(is.finite(fc_ardf$.estimate)))
  expect_true(all(fc_ardf$.error >= 0))
})

test_that("GPDF model works with basic functionality", {
  skip_on_cran()

  # Simple test data
  test_data <- data.frame(
    y = c(
      5,
      7,
      9,
      8,
      10,
      11,
      13,
      12,
      14,
      15,
      17,
      16,
      18,
      19,
      21,
      20,
      22,
      23,
      25,
      24
    ),
    season = rep(1:4, 5),
    time = 1:20
  )

  # Fit simple model
  mod <- SW(ffc_gam(
    y ~ fts(time, mean_only = TRUE, time_k = 6, k = 2) +
      fts(season, k = 3, time_k = 8),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  # Test GPDF forecasting
  newdata <- data.frame(
    season = c(1, 2),
    time = c(21, 22)
  )

  # Basic GPDF forecast
  fc_gpdf <- SW(forecast(
    mod,
    newdata = newdata,
    model = "GPDF",
    K = 1,
    chains = 1,
    iter = 300,
    summary = TRUE
  ))

  expect_true(inherits(fc_gpdf, "tbl_df"))
  expect_true(all(c(".estimate", ".error") %in% names(fc_gpdf)))
  expect_equal(nrow(fc_gpdf), nrow(newdata))
  expect_true(all(is.finite(fc_gpdf$.estimate)))
  expect_true(all(fc_gpdf$.error >= 0))
})

test_that("VARDF model works with basic functionality", {
  skip_on_cran()

  # Simple test data with more observations for VAR stability
  test_data <- data.frame(
    y = c(
      8,
      10,
      12,
      11,
      13,
      14,
      16,
      15,
      17,
      18,
      20,
      19,
      21,
      22,
      24,
      23,
      25,
      26,
      28,
      27,
      29,
      30,
      32,
      31,
      33,
      34,
      36,
      35,
      37,
      38
    ),
    season = rep(1:4, 9)[1:30],
    time = 1:30
  )

  # Fit model with multiple basis functions for VAR
  mod <- SW(ffc_gam(
    y ~ fts(time, mean_only = TRUE, time_k = 10, k = 3) +
      fts(season, k = 3, time_k = 8),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  # Test VARDF forecasting
  newdata <- data.frame(
    season = c(1, 2),
    time = c(31, 32)
  )

  # Basic VARDF forecast
  fc_vardf <- SW(forecast(
    mod,
    newdata = newdata,
    model = "VARDF",
    K = 2,
    lag = 1,
    chains = 1,
    iter = 300,
    summary = TRUE
  ))

  expect_true(inherits(fc_vardf, "tbl_df"))
  expect_true(all(c(".estimate", ".error") %in% names(fc_vardf)))
  expect_equal(nrow(fc_vardf), nrow(newdata))
  expect_true(all(is.finite(fc_vardf$.estimate)))
  expect_true(all(fc_vardf$.error >= 0))
})

test_that("Stan models validate parameters correctly", {
  skip_on_cran()

  # Simple test data
  test_data <- data.frame(
    y = c(
      5,
      7,
      9,
      8,
      10,
      11,
      13,
      12,
      14,
      15,
      17,
      16,
      18,
      19,
      21,
      20,
      22,
      23,
      25,
      24
    ),
    season = rep(1:4, 5),
    time = 1:20
  )

  mod <- SW(ffc_gam(
    y ~ fts(time, mean_only = TRUE, time_k = 5, k = 3) +
      fts(season, k = 3, time_k = 8),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  newdata <- data.frame(time = 21:22, season = 1:2)

  # Test invalid K parameter
  expect_error(
    forecast(mod, newdata = newdata, model = "ARDF", K = 0),
    "Assertion on 'K' failed: Must be >= 1"
  )

  # Test invalid lag parameter
  expect_error(
    forecast(mod, newdata = newdata, model = "ARDF", K = 1, lag = 0),
    "Assertion on 'lag' failed: Must be >= 1"
  )
})

test_that("ARDF forecasting works with mortality data example", {
  skip_on_cran()

  # Test exact example from ffc_gam documentation
  # This currently fails with duplicate .time column bug - needs to be fixed

  data("qld_mortality")

  # Use subset for faster testing
  qld_subset <- qld_mortality[qld_mortality$year >= 1990, ]

  mod <- SW(ffc_gam(
    deaths ~
      offset(log(population)) +
      sex +
      fts(age, k = 8, time_k = 10),
    time = "year",
    data = qld_subset,
    family = poisson(),
    engine = "bam"
  ))

  # Create forecast data exactly as in documentation example
  future_data <- expand.grid(
    age = unique(qld_subset$age),
    sex = unique(qld_subset$sex),
    year = 2021:2025,
    population = 1 # Use rate scale (deaths per person)
  )

  # This should work - currently fails due to duplicate .time column bug
  mortality_fc_ardf <- forecast(
    mod,
    newdata = future_data,
    model = "ARDF",
    K = 3,
    lag = 3,
    chains = 1,
    iter = 250,
    type = "expected"
  )

  # Test that results are reasonable
  expect_true(is.data.frame(mortality_fc_ardf))
  expect_true(nrow(mortality_fc_ardf) == nrow(future_data))
  expect_true(all(c(".estimate", ".error") %in% names(mortality_fc_ardf)))
  expect_true(all(mortality_fc_ardf$.estimate > 0)) # Mortality rates should be positive
})
