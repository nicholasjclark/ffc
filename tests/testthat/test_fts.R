# Tests for the fts() wrapper

test_that("fts() working correctly", {
  # Should return a list that correctly captures
  # the arguments, and that can create a valid mgcv smooth
  fts_obj <- fts(growth, k = 11, bs = "cr")
  expect_true(inherits(fts_obj, "list"))
  expect_true(fts_obj$term == "growth")
  expect_true(fts_obj$by == "NA")
  expect_true(fts_obj$time_bs == "ts")
  expect_true(fts_obj$time_k == 10)

  smooth_obj <- eval(parse(text = fts_obj$call))
  expect_true(inherits(smooth_obj, "cr.smooth.spec"))

  # Should error if using quoted variables
  expect_error(fts("growth", k = 11, bs = "cr"))

  # Also check that multidimensional bases are handled
  fts_obj <- fts(growth, age, k = 11, bs = "tp", time_bs = "ps")
  expect_true(inherits(fts_obj, "list"))
  expect_true(all.equal(fts_obj$term, c("growth", "age")))
  expect_true(fts_obj$by == "NA")
  expect_true(fts_obj$time_bs == "ps")
  expect_true(fts_obj$time_k == 10)

  smooth_obj <- eval(parse(text = fts_obj$call))
  expect_true(inherits(smooth_obj, "tensor.smooth.spec"))
  expect_true(all.equal(smooth_obj$term, c("growth", "age")))
})

test_that("ffc_gam handles response transformations correctly", {
  # Test various transformations of the response variable
  set.seed(4444)
  n <- 50

  # Create data with positive y for transformations
  test_data <- data.frame(
    y_raw = rgamma(n, shape = 2, rate = 1), # Positive values
    x = runif(n, 0, 1),
    time = 1:n
  )

  # Test sqrt transformation
  mod_sqrt <- SW(ffc_gam(
    sqrt(y_raw) ~ s(x, k = 6) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  expect_true(inherits(mod_sqrt, "ffc_gam"))
  expect_true(length(coef(mod_sqrt)) > 0)

  # Test log transformation
  mod_log <- SW(ffc_gam(
    log(y_raw) ~ s(x, k = 6) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  expect_true(inherits(mod_log, "ffc_gam"))
  expect_true(length(coef(mod_log)) > 0)

  # Test log2 transformation
  mod_log2 <- SW(ffc_gam(
    log2(y_raw) ~ s(x, k = 6) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  expect_true(inherits(mod_log2, "ffc_gam"))

  # Test log10 transformation
  mod_log10 <- SW(ffc_gam(
    log10(y_raw) ~ s(x, k = 6) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  expect_true(inherits(mod_log10, "ffc_gam"))

  # Test power transformation
  mod_power <- SW(ffc_gam(
    I(y_raw^2) ~ s(x, k = 6) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  expect_true(inherits(mod_power, "ffc_gam"))

  # Test predictions work with transformed responses
  newdata <- data.frame(
    y_raw = rgamma(5, shape = 2, rate = 1), # Not used in prediction
    x = runif(5, 0, 1),
    time = (n + 1):(n + 5)
  )

  # Predictions should work (though they're on the transformed scale)
  pred_sqrt <- predict(mod_sqrt, newdata = newdata)
  expect_equal(length(pred_sqrt), nrow(newdata))
  expect_true(all(is.finite(pred_sqrt)))

  pred_log <- predict(mod_log, newdata = newdata)
  expect_equal(length(pred_log), nrow(newdata))
  expect_true(all(is.finite(pred_log)))
})

test_that("ffc_gam handles response transformations with different families", {
  # Test response transformations with non-Gaussian families
  set.seed(5555)
  n <- 60

  # Create data suitable for different families
  test_data <- data.frame(
    count_raw = rpois(n, lambda = 5),
    gamma_raw = rgamma(n, shape = 3, rate = 1),
    x = runif(n, 0, 1),
    time = 1:n
  )

  # Test sqrt transformation with Poisson family (unusual but valid)
  # Note: this models sqrt(count) as if it were normally distributed
  mod_sqrt_pois <- SW(ffc_gam(
    sqrt(count_raw) ~ s(x, k = 6) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = gaussian() # Model transformed response as Gaussian
  ))

  expect_true(inherits(mod_sqrt_pois, "ffc_gam"))

  # Test log transformation with Gamma family
  # Note: this models log(gamma_var) as if it were normally distributed
  mod_log_gamma <- SW(ffc_gam(
    log(gamma_raw) ~ s(x, k = 6) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  expect_true(inherits(mod_log_gamma, "ffc_gam"))

  # Test with small values to avoid overflow issues
  test_data$small_vals <- runif(n, 0.1, 2)
  mod_exp <- SW(ffc_gam(
    exp(small_vals) ~ s(x, k = 6) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = Gamma() # Gamma family for positive continuous response
  ))

  expect_true(inherits(mod_exp, "ffc_gam"))

  # Test predictions with transformed response models
  newdata <- data.frame(
    count_raw = rpois(5, lambda = 5),
    gamma_raw = rgamma(5, shape = 3, rate = 1),
    small_vals = runif(5, 0.1, 2),
    x = runif(5, 0, 1),
    time = (n + 1):(n + 5)
  )

  pred_sqrt_pois <- predict(mod_sqrt_pois, newdata = newdata)
  expect_equal(length(pred_sqrt_pois), nrow(newdata))

  pred_log_gamma <- predict(mod_log_gamma, newdata = newdata)
  expect_equal(length(pred_log_gamma), nrow(newdata))
})

test_that("response transformations work with complex expressions", {
  # Test more complex response transformations
  set.seed(6666)
  n <- 40

  test_data <- data.frame(
    y1 = rnorm(n, 5, 1), # Mean around 5
    y2 = rnorm(n, 3, 0.5), # Mean around 3
    x = runif(n, 0, 1),
    time = 1:n
  )

  # Test arithmetic combinations
  mod_sum <- SW(ffc_gam(
    I(y1 + y2) ~ s(x, k = 6) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  expect_true(inherits(mod_sum, "ffc_gam"))

  # Test scaled response
  mod_scaled <- SW(ffc_gam(
    I(y1 * 2 + 1) ~ s(x, k = 6) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  expect_true(inherits(mod_scaled, "ffc_gam"))

  # Test with absolute values (to ensure positive)
  test_data$mixed_vals <- rnorm(n, 0, 2)
  mod_abs <- SW(ffc_gam(
    abs(mixed_vals) ~ s(x, k = 6) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = Gamma()
  ))

  expect_true(inherits(mod_abs, "ffc_gam"))

  # Verify predictions work
  newdata <- data.frame(
    y1 = rnorm(3, 5, 1),
    y2 = rnorm(3, 3, 0.5),
    mixed_vals = rnorm(3, 0, 2),
    x = runif(3, 0, 1),
    time = (n + 1):(n + 3)
  )

  pred_sum <- predict(mod_sum, newdata = newdata)
  expect_equal(length(pred_sum), nrow(newdata))

  pred_scaled <- predict(mod_scaled, newdata = newdata)
  expect_equal(length(pred_scaled), nrow(newdata))

  pred_abs <- predict(mod_abs, newdata = newdata)
  expect_equal(length(pred_abs), nrow(newdata))
})

test_that("response transformations work with arithmetic expressions inside functions", {
  # Test transformations with arithmetic expressions inside function calls
  # These were previously failing before delegation to mgcv
  set.seed(7777)
  n <- 45

  test_data <- data.frame(
    accel = rnorm(n, -50, 30), # Can be negative like motorcycle data
    pos_val = rgamma(n, 2, 1), # Always positive
    x = runif(n, 0, 1),
    time = 1:n
  )

  # Test log with offset (like log(accel + 150) from vignette)
  mod_log_offset <- SW(ffc_gam(
    log(accel + 200) ~ s(x, k = 6) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  expect_true(inherits(mod_log_offset, "ffc_gam"))
  expect_true(length(coef(mod_log_offset)) > 0)

  # Test sqrt with arithmetic expression
  mod_sqrt_expr <- SW(ffc_gam(
    sqrt(abs(accel) + 1) ~ s(x, k = 6) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  expect_true(inherits(mod_sqrt_expr, "ffc_gam"))
  expect_true(length(coef(mod_sqrt_expr)) > 0)

  # Test log10 with multiplication and addition
  mod_log10_complex <- SW(ffc_gam(
    log10(pos_val * 2 + 0.1) ~ s(x, k = 6) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = gaussian()
  ))

  expect_true(inherits(mod_log10_complex, "ffc_gam"))
  expect_true(length(coef(mod_log10_complex)) > 0)

  # Test distributional regression with complex response transformation
  mod_gaulss_transform <- SW(ffc_gam(
    list(
      log(accel + 200) ~ fts(time, k = 5, time_k = 5),
      ~ fts(time, k = 3, time_k = 5)
    ),
    family = mgcv::gaulss(),
    data = test_data,
    time = "time"
  ))

  expect_true(inherits(mod_gaulss_transform, "ffc_gam"))
  expect_true(inherits(mod_gaulss_transform, "ffc_gam_multi"))

  # Test predictions work with complex transformations
  newdata <- data.frame(
    accel = rnorm(4, -50, 30),
    pos_val = rgamma(4, 2, 1),
    x = runif(4, 0, 1),
    time = (n + 1):(n + 4)
  )

  pred_log_offset <- predict(mod_log_offset, newdata = newdata)
  pred_sqrt_expr <- predict(mod_sqrt_expr, newdata = newdata)
  pred_gaulss <- predict(mod_gaulss_transform, newdata = newdata)

  expect_equal(length(pred_log_offset), nrow(newdata))
  expect_equal(length(pred_sqrt_expr), nrow(newdata))
  expect_equal(nrow(pred_gaulss), nrow(newdata)) # Should return location parameter only

  # Ensure predictions are finite
  expect_true(all(is.finite(pred_log_offset)))
  expect_true(all(is.finite(pred_sqrt_expr)))
  expect_true(all(is.finite(pred_gaulss)))
})
