# Test post-processing of the example model

test_that("ffc_gam() is structured correctly", {
  expect_true(inherits(example_mod, "ffc_gam"))
  expect_true(inherits(example_mod, "bam"))
  expect_true(example_mod$time_var == "age_yr")
  eval(parse(text = example_mod$fts_smooths[[1]]$call))
})

test_that("mgcv post-processing works correctly", {
  expect_no_error(capture_output(summary(example_mod)))
  expect_no_error(capture_output(residuals(example_mod)))
  expect_no_error(capture_output(coef(example_mod)))
  expect_no_error(capture_output(predict(example_mod)))
  expect_no_error(capture_output(plot(example_mod, pages = 1)))
})

test_that("fts_coefs() works correctly", {
  functional_coefs <- fts_coefs(example_mod)
  expect_true(inherits(functional_coefs, "fts_ts"))
  expect_true(inherits(functional_coefs, "tbl_df"))

  # Times should match exactly to those included in the object
  expect_true(
    all.equal(
      sort(unique(functional_coefs$.time)),
      sort(unique(example_mod$model[[example_mod$time_var]]))
    )
  )

  expect_true(
    all.equal(
      sort(unique(functional_coefs[[example_mod$time_var]])),
      sort(unique(example_mod$model[[example_mod$time_var]]))
    )
  )

  # SEs should not be negative
  expect_range(functional_coefs$.se, lower = 0)
})
