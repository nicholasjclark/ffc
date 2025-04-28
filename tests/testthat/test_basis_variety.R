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
          time_bs = 'cr',
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
        bs = "mrf", xt = xt,
        k = 6,
        time_k = 3
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
