# Test post-processing of the example model

test_that("ffc_gam() is structured correctly", {
  expect_true(inherits(example_mod, "ffc_gam"))
  expect_true(inherits(example_mod, "bam"))
  expect_true(example_mod$time_var == "age_yr")
  eval(parse(text = example_mod$fts_smooths[[1]]$call))

  expect_true(
    identical(
      attr(example_mod2$gam_init[[1]], "knots"),
      list(season = c(0.5, 12.5))
    )
  )
})

test_that("mgcv post-processing works correctly", {
  expect_no_error(capture_output(summary(example_mod)))
  expect_no_error(capture_output(residuals(example_mod)))
  expect_no_error(capture_output(coef(example_mod)))
  expect_no_error(capture_output(predict(example_mod)))
  expect_no_error(capture_output(plot(example_mod, pages = 1)))

  expect_no_error(capture_output(summary(example_mod2)))
  expect_no_error(capture_output(residuals(example_mod2)))
  expect_no_error(capture_output(coef(example_mod2)))
  expect_no_error(capture_output(predict(example_mod2)))
  expect_no_error(capture_output(plot(example_mod2, pages = 1)))
})

test_that("fts_coefs() works correctly", {
  functional_coefs <- fts_coefs(example_mod)
  expect_true(inherits(functional_coefs, "fts_ts"))
  expect_true(inherits(functional_coefs, "tbl_df"))
  expect_true(
    all.equal(
      attr(functional_coefs, "time_var"),
      example_mod$time_var
    )
  )
  expect_true(
    attr(functional_coefs, "summarized")
  )

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

  # Plot should return a ggplot
  expect_ggplot(
    autoplot(functional_coefs)
  )

  # SEs should not be negative
  expect_range(functional_coefs$.se, lower = 0)

  # Also need to test outputs when summary = FALSE
  functional_coefs <- fts_coefs(
    example_mod,
    summary = FALSE,
    times = 3
  )
  expect_true(inherits(functional_coefs, "fts_ts"))
  expect_true(inherits(functional_coefs, "tbl_df"))
  expect_true(
    all.equal(
      attr(functional_coefs, "time_var"),
      example_mod$time_var
    )
  )
  expect_false(
    attr(functional_coefs, "summarized")
  )
  expect_true(
    length(unique(functional_coefs$.realisation)) == 3
  )
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
  expect_ggplot(
    autoplot(functional_coefs)
  )

  # And ensure only a single fts basis from the example_mod2
  functional_coefs <- fts_coefs(
    example_mod2,
    summary = FALSE,
    times = 3
  )
  expect_true(
    length(unique(functional_coefs$.basis)) == 1L
  )
  expect_ggplot(
    autoplot(functional_coefs)
  )
})
