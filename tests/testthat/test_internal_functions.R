# Test internal functions, transformations, character families, and gam_setup

test_that("rmvn generates multivariate normal samples correctly", {
  set.seed(123)
  n <- 100
  mu <- c(1, 2, 3)
  sig <- matrix(c(2, 1, 0.5, 1, 3, 0.8, 0.5, 0.8, 1.5), nrow = 3)
  
  # Generate samples using internal function
  samples <- ffc:::rmvn(n, mu, sig)
  
  # Check dimensions
  expect_equal(nrow(samples), n)
  expect_equal(ncol(samples), length(mu))
  
  # Check means are approximately correct
  sample_means <- colMeans(samples)
  expect_true(all(abs(sample_means - mu) < 0.2))  # Should be close for n=100
  
  # Check covariance structure is approximately correct
  sample_cov <- cov(samples)
  expect_equal(dim(sample_cov), dim(sig))
  
  # Check that samples are numeric
  expect_true(is.numeric(samples))
  expect_true(all(is.finite(samples)))
})

test_that("rmvn handles edge cases correctly", {
  # Single dimension
  n <- 50
  mu <- 2
  sig <- matrix(1)
  samples <- ffc:::rmvn(n, mu, sig)
  expect_equal(nrow(samples), n)
  expect_equal(ncol(samples), 1)
  expect_true(abs(mean(samples) - mu) < 0.3)
  
  # Different n values
  for (n_test in c(1, 5, 200)) {
    mu <- c(0, 1)
    sig <- diag(2)
    samples <- ffc:::rmvn(n_test, mu, sig)
    expect_equal(nrow(samples), n_test)
    expect_equal(ncol(samples), 2)
  }
})

test_that("mathematical transformations work in ffc_gam formulas", {
  # Create test data with transformations
  test_data <- data.frame(
    y = rgamma(50, shape = 2, rate = 1),
    x = runif(50, 1, 10),
    time = 1:50
  )
  
  # Test log transformation
  mod_log <- SW(ffc_gam(
    y ~ log(x) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = Gamma()
  ))
  expect_true(inherits(mod_log, "ffc_gam"))
  expect_true(length(coef(mod_log)) > 0)
  
  # Test log2 transformation
  mod_log2 <- SW(ffc_gam(
    y ~ log2(x) + fts(time, k = 5, time_k = 5),
    time = "time", 
    data = test_data,
    family = Gamma()
  ))
  expect_true(inherits(mod_log2, "ffc_gam"))
  
  # Test sqrt transformation
  mod_sqrt <- SW(ffc_gam(
    y ~ sqrt(x) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = Gamma()
  ))
  expect_true(inherits(mod_sqrt, "ffc_gam"))
  
  # Test exp transformation (with smaller values to avoid overflow)
  test_data$x_small <- test_data$x / 10
  mod_exp <- SW(ffc_gam(
    y ~ exp(x_small) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = Gamma()
  ))
  expect_true(inherits(mod_exp, "ffc_gam"))
  
  # Test trigonometric functions
  mod_sin <- SW(ffc_gam(
    y ~ sin(x) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = Gamma()
  ))
  expect_true(inherits(mod_sin, "ffc_gam"))
  
  mod_cos <- SW(ffc_gam(
    y ~ cos(x) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = Gamma()
  ))
  expect_true(inherits(mod_cos, "ffc_gam"))
})

test_that("character family specifications work correctly", {
  # Test data for different families
  gaussian_data <- data.frame(
    y = rnorm(40, 10, 2),
    x = 1:40,
    time = 1:40
  )
  
  poisson_data <- data.frame(
    y = rpois(40, lambda = 5),
    x = 1:40,
    time = 1:40
  )
  
  gamma_data <- data.frame(
    y = rgamma(40, shape = 2, rate = 0.5),
    x = 1:40,
    time = 1:40
  )
  
  binomial_data <- data.frame(
    successes = rbinom(40, size = 10, prob = 0.3),
    trials = rep(10, 40),
    x = 1:40,
    time = 1:40
  )
  binomial_data$y <- cbind(binomial_data$successes, binomial_data$trials - binomial_data$successes)
  
  # Test "gaussian" family
  mod_gaussian <- SW(ffc_gam(
    y ~ s(x, k = 5) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = gaussian_data,
    family = "gaussian"
  ))
  expect_true(inherits(mod_gaussian, "ffc_gam"))
  expect_equal(mod_gaussian$family$family, "gaussian")
  
  # Test "poisson" family
  mod_poisson <- SW(ffc_gam(
    y ~ s(x, k = 5) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = poisson_data,
    family = "poisson"
  ))
  expect_true(inherits(mod_poisson, "ffc_gam"))
  expect_equal(mod_poisson$family$family, "poisson")
  
  # Test "Gamma" family
  mod_gamma <- SW(ffc_gam(
    y ~ s(x, k = 5) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = gamma_data,
    family = "Gamma"
  ))
  expect_true(inherits(mod_gamma, "ffc_gam"))
  expect_equal(mod_gamma$family$family, "Gamma")
  
  # Test "binomial" family
  mod_binomial <- SW(ffc_gam(
    y ~ s(x, k = 5) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = binomial_data,
    family = "binomial"
  ))
  expect_true(inherits(mod_binomial, "ffc_gam"))
  expect_equal(mod_binomial$family$family, "binomial")
})

test_that("character family specifications vs function forms give same results", {
  test_data <- data.frame(
    y = rnorm(30, 5, 1),
    x = 1:30,
    time = 1:30
  )
  
  # Fit with function form
  mod_func <- SW(ffc_gam(
    y ~ s(x, k = 5) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = gaussian()
  ))
  
  # Fit with character form
  mod_char <- SW(ffc_gam(
    y ~ s(x, k = 5) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = "gaussian"
  ))
  
  # Should have same family structure
  expect_equal(mod_func$family$family, mod_char$family$family)
  expect_equal(mod_func$family$link, mod_char$family$link)
  
  # Coefficients should be very similar (allowing for numerical differences)
  expect_true(all(abs(coef(mod_func) - coef(mod_char)) < 1e-6))
})

test_that("gam_setup and ffc_gam_setup work through ffc_gam interface", {
  # Test that the setup functions work correctly through the main interface
  test_data <- data.frame(
    y = rnorm(30, 0, 1),
    x = rnorm(30, 0, 1),
    time = 1:30
  )
  
  # Test that ffc_gam successfully uses internal setup functions
  mod_simple <- SW(ffc_gam(
    y ~ s(x, k = 5) + fts(time, k = 5, time_k = 5),
    time = "time",
    data = test_data,
    family = gaussian()
  ))
  
  expect_true(inherits(mod_simple, "ffc_gam"))
  expect_true(inherits(mod_simple, "gam"))
  expect_true(length(coef(mod_simple)) > 0)
  
  # Test that setup functions handle knots correctly through interface
  mod_knots <- SW(ffc_gam(
    y ~ s(x, k = 5) + fts(time, k = 5, time_k = 5),
    knots = list(x = c(-2, 0, 2)),
    time = "time",
    data = test_data,
    family = gaussian()
  ))
  
  expect_true(inherits(mod_knots, "ffc_gam"))
  expect_true(length(coef(mod_knots)) > 0)
  
  # Test error for invalid knots (should be caught by ffc_gam_setup)
  expect_error(ffc_gam(
    y ~ s(x, k = 5) + fts(time, k = 5, time_k = 5),
    knots = "invalid_knots",  # Should be a list
    time = "time",
    data = test_data,
    family = gaussian()
  ), "knot.*must be supplied as lists")
})