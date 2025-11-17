# Test plotting functionality for autoplot.fts_ts

test_that("autoplot works with single basis function (n_basis = 1L)", {
  # Use example_mod which has only one basis function
  functional_coefs <- fts_coefs(example_mod, summary = TRUE)
  
  # Check number of basis functions
  n_basis <- length(unique(functional_coefs$.basis))
  expect_equal(n_basis, 1L)
  
  # Test plot creation
  p <- autoplot(functional_coefs)
  expect_ggplot(p)
  
  # Check that it has a title (for single basis)
  expect_true(any(grepl("basis", p$labels$title, ignore.case = TRUE)))
  
  # Should have FacetNull for single basis (not faceted)
  expect_true(inherits(p$facet, "FacetNull"))
})

test_that("autoplot works with multiple basis functions (n_basis > 1L)", {
  # Create a model with multiple basis functions
  multi_basis_data <- data.frame(
    y = rnorm(100, 10, 2),
    x1 = runif(100, 0, 10),
    x2 = runif(100, -5, 5),
    time = 1:100
  )
  
  # Fit model with multiple fts terms to get multiple basis functions
  multi_mod <- SW(ffc_gam(
    y ~ fts(x1, k = 3, time_k = 5) + fts(x2, k = 4, time_k = 5),
    time = "time",
    data = multi_basis_data,
    family = gaussian()
  ))
  
  # Extract coefficients
  functional_coefs <- fts_coefs(multi_mod, summary = TRUE)
  
  # Check that we have multiple basis functions
  n_basis <- length(unique(functional_coefs$.basis))
  expect_true(n_basis > 1L)
  
  # Test plot creation
  p <- autoplot(functional_coefs)
  expect_ggplot(p)
  
  # Should have facets for multiple basis functions
  expect_true(inherits(p$facet, "FacetWrap"))
  expect_true("~.basis" %in% as.character(p$facet$params$facets))
  
  # Check facet parameters (scales parameter may be named differently)
  facet_params <- p$facet$params
  expect_true("scales" %in% names(facet_params) || "free" %in% names(facet_params))
  expect_equal(p$facet$params$ncol, min(4, n_basis))
  
  # Should not have a title (titles are only for single basis)
  expect_null(p$labels$title)
})

test_that("autoplot works with non-summarized coefficients (multiple realizations)", {
  # Test with non-summarized data (should have .realisation grouping)
  functional_coefs <- fts_coefs(example_mod, summary = FALSE, times = 3)
  
  # Check that it's not summarized
  expect_false(attr(functional_coefs, "summarized"))
  expect_true(".realisation" %in% names(functional_coefs))
  expect_true(length(unique(functional_coefs$.realisation)) == 3)
  
  # Test plot creation
  p <- autoplot(functional_coefs)
  expect_ggplot(p)
  
  # Should have grouping aesthetic for .realisation
  expect_true("group" %in% names(p$mapping))
  
  # Check that colour and group aesthetics are properly set
  expect_true(rlang::as_label(p$mapping$colour) == ".basis")
  expect_true(rlang::as_label(p$mapping$group) == ".realisation")
})

test_that("autoplot handles basis name formatting correctly", {
  # Create model and extract coefficients
  functional_coefs <- fts_coefs(example_mod, summary = TRUE)
  
  # Test plot creation
  p <- autoplot(functional_coefs)
  
  # Check that basis names are properly formatted
  plot_data <- p$data
  expect_true(all(grepl("basis", plot_data$.basis)))
  expect_false(any(grepl("fts_bs_", plot_data$.basis)))
  expect_false(any(grepl("fts_", plot_data$.basis, fixed = TRUE)))
})

test_that("autoplot works with different numbers of basis functions", {
  # Test various numbers of basis functions to check ncol logic
  
  # Test with 2 basis functions
  test_data_2 <- data.frame(
    y = rnorm(60, 5, 1),
    x1 = runif(60, 0, 5),
    x2 = runif(60, -3, 3),
    time = 1:60
  )
  
  mod_2basis <- SW(ffc_gam(
    y ~ fts(x1, k = 2, time_k = 4) + fts(x2, k = 3, time_k = 4),
    time = "time",
    data = test_data_2,
    family = gaussian()
  ))
  
  coefs_2 <- fts_coefs(mod_2basis, summary = TRUE)
  n_basis_2 <- length(unique(coefs_2$.basis))
  
  p2 <- autoplot(coefs_2)
  expect_ggplot(p2)
  expect_true(inherits(p2$facet, "FacetWrap"))
  expect_equal(p2$facet$params$ncol, min(4, n_basis_2))
  
  # Test with more basis functions (should cap at 4 columns)
  test_data_5 <- data.frame(
    y = rnorm(100, 3, 0.5),
    x1 = runif(100, 0, 2),
    x2 = runif(100, -1, 1),
    x3 = runif(100, 1, 3),
    time = 1:100
  )
  
  mod_5basis <- SW(ffc_gam(
    y ~ fts(x1, k = 3, time_k = 4) + 
        fts(x2, k = 2, time_k = 4) + 
        fts(x3, k = 3, time_k = 4),
    time = "time", 
    data = test_data_5,
    family = gaussian()
  ))
  
  coefs_5 <- fts_coefs(mod_5basis, summary = TRUE)
  n_basis_5 <- length(unique(coefs_5$.basis))
  
  p5 <- autoplot(coefs_5)
  expect_ggplot(p5)
  expect_true(inherits(p5$facet, "FacetWrap"))
  # Should cap at 4 columns regardless of number of basis functions
  expect_equal(p5$facet$params$ncol, min(4, n_basis_5))
})

test_that("autoplot aesthetic mappings are correct", {
  # Test with summarized data
  coefs_summ <- fts_coefs(example_mod, summary = TRUE)
  p_summ <- autoplot(coefs_summ)
  
  # Should have x, y, colour but not group for summarized data
  expect_true("x" %in% names(p_summ$mapping))
  expect_true("y" %in% names(p_summ$mapping))
  expect_true("colour" %in% names(p_summ$mapping))
  expect_false("group" %in% names(p_summ$mapping))
  
  # Check specific mappings
  expect_true(rlang::as_label(p_summ$mapping$x) == attr(coefs_summ, "time_var"))
  expect_true(rlang::as_label(p_summ$mapping$y) == ".estimate")
  expect_true(rlang::as_label(p_summ$mapping$colour) == ".basis")
  
  # Test with non-summarized data
  coefs_raw <- fts_coefs(example_mod, summary = FALSE, times = 2)
  p_raw <- autoplot(coefs_raw)
  
  # Should have x, y, colour AND group for non-summarized data
  expect_true("x" %in% names(p_raw$mapping))
  expect_true("y" %in% names(p_raw$mapping))
  expect_true("colour" %in% names(p_raw$mapping))
  expect_true("group" %in% names(p_raw$mapping))
  
  # Check group mapping for non-summarized
  expect_true(rlang::as_label(p_raw$mapping$group) == ".realisation")
})

test_that("autoplot handles geom_line vs geom_borderline correctly", {
  # This tests the conditional logic for ggborderline package
  functional_coefs <- fts_coefs(example_mod, summary = TRUE)
  
  # Test plot creation (should work regardless of ggborderline availability)
  expect_no_error({
    p <- autoplot(functional_coefs)
  })
  expect_ggplot(p)
  
  # Should have at least one layer (either geom_line or geom_borderline)
  expect_true(length(p$layers) >= 1)
  
  # Should have viridis color scale
  scale_names <- sapply(p$scales$scales, function(x) class(x)[1])
  expect_true(any(grepl("ScaleDiscrete", scale_names)))
})

test_that("autoplot labeller works with long basis names", {
  # Create model with potentially long basis names
  test_data <- data.frame(
    y = rnorm(60, 2, 0.3),
    very_long_variable_name = runif(60, 0, 10),
    another_long_name = runif(60, -2, 2),
    time = 1:60
  )
  
  mod_long <- SW(ffc_gam(
    y ~ fts(very_long_variable_name, k = 2, time_k = 4) + 
        fts(another_long_name, k = 2, time_k = 4),
    time = "time",
    data = test_data,
    family = gaussian()
  ))
  
  coefs_long <- fts_coefs(mod_long, summary = TRUE)
  p_long <- autoplot(coefs_long)
  
  expect_ggplot(p_long)
  expect_true(inherits(p_long$facet, "FacetWrap"))
  
  # Check that labeller is set up (structure may vary)
  expect_true(!is.null(p_long$facet$params$labeller))
})