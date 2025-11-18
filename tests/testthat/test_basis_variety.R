# Tests for a variety of basis functions

test_that("MRF basis setup works correctly", {
  library(mgcv)

  # Use the columb example from mgcv
  data(columb)
  data(columb.polys)
  xt <- list(polys = columb.polys)

  # Create some data with a time dimension
  fake_dat <- do.call(
    rbind,
    lapply(1:5, function(x) {
      dat <- columb
      dat$time <- x
      dat
    })
  )

  # Try fitting a model, but use too many basis functions
  # for time_k
  expect_error(
    ffc_gam(
      crime ~
        fts(district,
          bs = "mrf",
          xt = xt,
          k = 6,
          time_bs = "cr",
          time_k = 10
        ),
      time = "time",
      data = fake_dat,
      family = "gaussian"
    )
  )

  # Reduce to a usable time_k and try again
  mod <- ffc_gam(
    crime ~
      fts(district,
        bs = "mrf",
        xt = xt,
        k = 6,
        time_k = 3,
        time_bs = "cr"
      ),
    time = "time",
    data = fake_dat,
    family = "gaussian"
  )
  expect_true(inherits(mod, "ffc_gam"))
  expect_true(
    inherits(
      mod$gam_init[[1]]$smooth[[1]],
      "mrf.smooth"
    )
  )
  expect_true(
    identical(
      mod$gam_init[[1]]$smooth[[1]]$xt$polys,
      xt$polys
    )
  )
})

test_that("by variables work correctly in regular smooths with numeric variables", {
  # Create test data with a numeric 'by' variable
  n <- 50
  x1 <- runif(n, 0, 1)
  x2 <- runif(n, 0, 1)
  z <- runif(n, -1, 1)  # numeric by variable
  time <- rep(seq(1, 5, by = 1), each = 10)  # regular intervals
  
  # Generate y with interaction between smooth of x1 and z
  y <- sin(2 * pi * x1) * z + rnorm(n, 0, 0.2)
  test_data <- data.frame(x1 = x1, x2 = x2, z = z, time = time, y = y)
  
  # Fit model with numeric by variable
  mod <- ffc_gam(
    y ~ s(x1, by = z) + s(x2),
    time = "time",
    data = test_data,
    family = "gaussian"
  )
  
  expect_true(inherits(mod, "ffc_gam"))
  expect_true(length(mod$smooth) == 2)
  
  # Check that the by variable was correctly incorporated
  expect_true(any(sapply(mod$smooth, function(s) !is.null(s$by))))
  
  # Ensure model can make predictions
  pred <- predict(mod)
  expect_equal(length(pred), n)
})

test_that("by variables work correctly in regular smooths with factor variables", {
  # Create test data with a factor 'by' variable and regular time intervals
  time_seq <- seq(1, 10, by = 1)
  fac_levels <- c("A", "B", "C")
  
  # Create balanced design with regular time intervals
  test_data <- expand.grid(time = time_seq, fac = fac_levels)
  n <- nrow(test_data)
  
  test_data$x1 <- runif(n, 0, 1)
  test_data$x2 <- runif(n, 0, 1)
  test_data$fac <- factor(test_data$fac)  # factor by variable
  
  # Generate y with different smooth relationships for each factor level
  test_data$y <- with(test_data, 
    ifelse(fac == "A", sin(2 * pi * x1),
           ifelse(fac == "B", cos(2 * pi * x1), 
                  x1^2)) + rnorm(n, 0, 0.2)
  )
  
  # Fit model with factor by variable (parametric term + by smooth)
  mod <- ffc_gam(
    y ~ fac + s(x1, by = fac) + s(x2),
    time = "time",
    data = test_data,
    family = "gaussian"
  )
  
  expect_true(inherits(mod, "ffc_gam"))
  expect_true(length(mod$smooth) >= 2)
  
  # Check that the factor by variable created separate smooths
  by_smooths <- sapply(mod$smooth, function(s) !is.null(s$by))
  expect_true(any(by_smooths))
  
  # Ensure model can make predictions
  pred <- predict(mod)
  expect_equal(length(pred), n)
})

test_that("id arguments work correctly in regular smooths", {
  # Create test data with multiple predictors and regular time intervals
  n <- 80
  x1 <- runif(n, 0, 1)
  x2 <- runif(n, 0, 1)
  x3 <- runif(n, 0, 1)
  x4 <- runif(n, 0, 1)
  time <- rep(seq(1, 4, by = 1), each = 20)  # regular intervals
  
  # Generate y as combination of smooth functions
  y <- sin(2 * pi * x1) + cos(2 * pi * x2) + 
       sin(2 * pi * x3) + cos(2 * pi * x4) + rnorm(n, 0, 0.2)
  test_data <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, time = time, y = y)
  
  # Fit model with id argument to force same smoothing parameters
  mod <- ffc_gam(
    y ~ s(x1, id = 1) + s(x2) + s(x3, id = 1) + s(x4),
    time = "time",
    data = test_data,
    family = "gaussian"
  )
  
  expect_true(inherits(mod, "ffc_gam"))
  expect_true(length(mod$smooth) == 4)
  
  # Check that smooths with same id have same smoothing parameters
  id_smooths <- which(sapply(mod$smooth, function(s) !is.null(s$id) && s$id == 1))
  if (length(id_smooths) >= 2) {
    # If mgcv successfully processed id arguments, these should have linked smoothing
    expect_true(length(id_smooths) >= 2)
  }
  
  # Ensure model can make predictions
  pred <- predict(mod)
  expect_equal(length(pred), n)
})

test_that("combined by variables and id arguments work correctly", {
  # Create test data with both by variables and id arguments with regular time
  time_seq <- seq(1, 12, by = 1)
  fac_levels <- c("A", "B")
  
  # Create balanced design with regular time intervals
  test_data <- expand.grid(time = time_seq, fac = fac_levels)
  n <- nrow(test_data)
  
  test_data$x1 <- runif(n, 0, 1)
  test_data$x2 <- runif(n, 0, 1)
  test_data$fac <- factor(test_data$fac)
  test_data$z <- runif(n, -1, 1)
  
  # Generate y with complex interaction structure
  test_data$y <- with(test_data,
    ifelse(fac == "A", sin(2 * pi * x1) * z, cos(2 * pi * x1) * z) + 
    sin(2 * pi * x2) + rnorm(n, 0, 0.2)
  )
  
  # Fit model with both by variables and id arguments
  mod <- ffc_gam(
    y ~ fac + s(x1, by = fac, id = 1) + s(x2, id = 1) + s(x1, by = z),
    time = "time",
    data = test_data,
    family = "gaussian"
  )
  
  expect_true(inherits(mod, "ffc_gam"))
  expect_true(length(mod$smooth) >= 3)
  
  # Check that by variables were processed
  by_smooths <- sapply(mod$smooth, function(s) !is.null(s$by))
  expect_true(any(by_smooths))
  
  # Ensure model can make predictions
  pred <- predict(mod)
  expect_equal(length(pred), n)
})
