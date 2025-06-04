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
  fts_obj <- fts(growth, age,
    k = 11, bs = "tp",
    time_bs = "ps"
  )
  expect_true(inherits(fts_obj, "list"))
  expect_true(all.equal(fts_obj$term, c("growth", "age")))
  expect_true(fts_obj$by == "NA")
  expect_true(fts_obj$time_bs == "ps")
  expect_true(fts_obj$time_k == 10)

  smooth_obj <- eval(parse(text = fts_obj$call))
  expect_true(inherits(smooth_obj, "tensor.smooth.spec"))
  expect_true(all.equal(smooth_obj$term, c("growth", "age")))
})
