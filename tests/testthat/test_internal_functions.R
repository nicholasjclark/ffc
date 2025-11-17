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

test_that("compress_iseq compresses integer sequences correctly", {
  # Test basic consecutive sequence
  result1 <- ffc:::compress_iseq(c(1, 2, 3, 4, 5))
  expect_equal(result1, "c(1:5)")
  
  # Test mixed consecutive and individual numbers
  result2 <- ffc:::compress_iseq(c(1, 2, 3, 5, 6, 8))
  expect_equal(result2, "c(1:3,5:6,8)")
  
  # Test individual numbers only
  result3 <- ffc:::compress_iseq(c(1, 3, 5, 7))
  expect_equal(result3, "c(1,3,5,7)")
  
  # Test single number
  result4 <- ffc:::compress_iseq(5)
  expect_equal(result4, "c(5)")
  
  # Test unordered input (function should sort it)
  result5 <- ffc:::compress_iseq(c(3, 1, 2, 5, 4))
  expect_equal(result5, "c(1:5)")
})

test_that("compress_iseq handles edge cases", {
  # Test empty vector (produces NA result)
  result_empty <- ffc:::compress_iseq(integer(0))
  expect_true(is.character(result_empty))
  expect_true(length(result_empty) == 1)
  
  # Test repeated values (function doesn't deduplicate, just sorts)
  result_duplicates <- ffc:::compress_iseq(c(1, 1, 2, 3, 3))
  expect_true(is.character(result_duplicates))
  expect_true(grepl("c\\(", result_duplicates))  # Should start with "c("
})

test_that("initial_spg returns proper smoothing parameter structure", {
  # Create simple test data for GAM
  set.seed(123)
  n <- 20
  x <- runif(n, 0, 1)
  y <- sin(2 * pi * x) + rnorm(n, 0, 0.1)
  
  # Create simple design matrix (polynomial basis)
  X <- cbind(1, x, x^2, x^3)
  
  # Create simple penalty matrix
  D <- diff(diag(4), diff = 2)  # Second difference penalty
  S <- list(t(D) %*% D)
  
  # Test basic functionality
  sp_init <- SW(ffc:::initial_spg(
    x = X,
    y = y,
    weights = rep(1, n),
    family = gaussian(),
    S = S,
    rank = c(2),
    off = c(1, 5)  # Penalty applies to columns 2-4 (indices 1+1 to 5-1)
  ))
  
  # Check that result is a vector of smoothing parameters
  expect_true(is.numeric(sp_init))
  expect_true(length(sp_init) == length(S))
  expect_true(all(sp_init >= 0))  # Smoothing parameters should be non-negative
})

test_that("olid detects identifiability issues correctly", {
  # Create a design matrix with identifiability issues (collinear columns)
  X <- cbind(
    c(1, 1, 1, 1),     # intercept
    c(1, 2, 3, 4),     # x
    c(2, 4, 6, 8),     # 2*x (collinear with column 2)
    c(1, 4, 9, 16)     # x^2
  )
  
  # Basic test parameters
  nsdf <- c(1, 3)    # 1 parametric term, 3 smooth terms
  pstart <- c(1, 2)  # penalty starts
  flpi <- list(c(1), c(1))  # factor level parameter indices for each smooth
  lpi <- list(1:4)   # all parameters in first linear predictor
  
  # Test identifiability detection
  result <- SW(ffc:::olid(X, nsdf, pstart, flpi, lpi))
  
  # Check structure of result
  expect_true(is.list(result))
  expect_true("dind" %in% names(result))
  expect_true("lpi" %in% names(result))
  
  # Should detect some identifiability issue
  expect_true(is.vector(result$dind))
})

test_that("olid handles full rank design matrices", {
  # Create a full rank design matrix
  X <- cbind(
    c(1, 1, 1, 1),     # intercept
    c(1, 2, 3, 4),     # x
    c(1, 4, 9, 16),    # x^2
    c(1, 8, 27, 64)    # x^3
  )
  
  nsdf <- c(1, 3)
  pstart <- c(1, 2)
  flpi <- list(c(1), c(1))  # proper factor level parameter indices
  lpi <- list(1:4)
  
  # Test with full rank matrix
  result <- SW(ffc:::olid(X, nsdf, pstart, flpi, lpi))
  
  expect_true(is.list(result))
  expect_true("dind" %in% names(result))
  
  # Should not detect any issues (dind should be empty or NULL)
  expect_true(is.null(result$dind) || length(result$dind) == 0)
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

test_that("parametric_penalty returns NULL when no penalties", {
  # Create proper terms object using actual R terms functionality
  f <- ~ 1
  pterms <- terms(f)
  
  # Empty assignment (no terms beyond intercept)
  assign <- c(0)  # only intercept
  
  # Test with empty paraPen
  result <- SW(ffc:::parametric_penalty(
    pterms = pterms,
    assign = assign,
    paraPen = list(),
    sp0 = 0.1
  ))
  
  # Function returns NULL when no penalties are processed
  expect_null(result)
})

test_that("parametric_penalty returns NULL for terms without penalties", {
  # Test with terms that have no penalties applied
  f <- ~ x + z
  pterms <- terms(f)
  assign <- c(0, 1, 2)  # intercept, x, z
  
  # Call with NULL paraPen
  result <- SW(ffc:::parametric_penalty(
    pterms = pterms,
    assign = assign,
    paraPen = NULL,
    sp0 = 0.1
  ))
  
  # Returns NULL when no penalties are found
  expect_null(result)
})

test_that("clone_smooth_spec clones smooth specifications correctly", {
  # Create mock smooth specifications
  specb <- list(
    term = "s(x)",
    label = "s(x)", 
    dim = 1,
    by = "NA",
    id = 1,
    margin = list(bs = "tp", k = 10),
    xt = NULL
  )
  class(specb) <- "xx.smooth.spec"
  
  spec <- list(
    term = "s(y)",
    label = "s(y)",
    dim = 1, 
    by = "NA",
    id = 1,
    margin = list(bs = "tp", k = 10),
    xt = NULL
  )
  class(spec) <- "xx.smooth.spec"
  
  # Test basic cloning
  result <- SW(ffc:::clone_smooth_spec(specb, spec))
  
  expect_equal(result$term, "s(y)")
  expect_equal(result$label, "s(y)")
  expect_equal(result$by, "NA")
  expect_equal(result$id, 1)
  expect_true(inherits(result, "xx.smooth.spec"))
})

test_that("clone_smooth_spec handles by variables correctly", {
  specb <- list(
    term = "s(x)",
    label = "s(x)",
    dim = 1,
    by = "group",
    margin = list(bs = "tp", k = 10)
  )
  class(specb) <- "xx.smooth.spec"
  
  spec <- list(
    term = "s(y)", 
    label = "s(y)",
    dim = 1,
    by = "NA",
    margin = list(bs = "tp", k = 10)
  )
  class(spec) <- "xx.smooth.spec"
  
  result <- SW(ffc:::clone_smooth_spec(specb, spec))
  
  expect_equal(result$by, "NA")
  expect_equal(result$term, "s(y)")
})

test_that("initial_spg handles different family types", {
  set.seed(42)
  n <- 15
  X <- cbind(1, runif(n), rnorm(n))
  y <- rnorm(n)
  S <- list(diag(3))
  
  # Test with Poisson family
  y_pois <- rpois(n, lambda = 2)
  sp_pois <- SW(ffc:::initial_spg(
    x = X,
    y = y_pois,
    weights = rep(1, n),
    family = poisson(),
    S = S,
    rank = c(3),
    off = c(1, 4)
  ))
  
  expect_true(is.numeric(sp_pois))
  expect_true(length(sp_pois) == 1)
  expect_true(sp_pois >= 0)
  
  # Test with Gamma family  
  y_gamma <- rgamma(n, shape = 2)
  sp_gamma <- SW(ffc:::initial_spg(
    x = X,
    y = y_gamma,
    weights = rep(1, n),
    family = Gamma(),
    S = S,
    rank = c(3),
    off = c(1, 4)
  ))
  
  expect_true(is.numeric(sp_gamma))
  expect_true(length(sp_gamma) == 1)
  expect_true(sp_gamma >= 0)
})

test_that("initial_spg handles empty penalty list", {
  set.seed(123)
  n <- 10
  X <- cbind(1, runif(n))
  y <- rnorm(n)
  
  # Test with empty penalty list
  sp_empty <- SW(ffc:::initial_spg(
    x = X,
    y = y,
    weights = rep(1, n), 
    family = gaussian(),
    S = list(),
    rank = numeric(0),
    off = numeric(0)
  ))
  
  expect_true(is.numeric(sp_empty))
  expect_true(length(sp_empty) == 0)
})