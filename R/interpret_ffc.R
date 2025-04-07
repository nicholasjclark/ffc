#' Interpret the formula specified to ffc and replace any fts() terms
#' with the correct time-varying coefficient terms
#' @importFrom stats formula terms as.formula terms.formula
#' @noRd
interpret_ffc <- function(
    formula,
    data,
    time_var = "time",
    gam_init = list(),
    newdata = NULL) {
  # Factors in formula
  facs <- colnames(attr(terms.formula(formula), "factors"))

  # Check if formula has an intercept
  keep_intercept <- attr(terms(formula), "intercept") == 1

  # Check if any terms use the fts() wrapper
  response <- terms.formula(formula)[[2]]
  tf <- attr(terms.formula(formula, keep.order = TRUE), "term.labels")
  which_dynamics <- grep("fts(", tf, fixed = TRUE)

  # Update the formula to the correct time-varying coefficient implementation
  fts_smooths <- list()
  if (length(which_dynamics) != 0L) {
    termlabs <- attr(terms(formula, keep.order = TRUE), "term.labels")
    fts_smooths <- vector(length = length(which_dynamics), mode = "list")

    # Evaluate the fts() terms and update the formula appropriately
    if (length(which_dynamics > 1)) {
      termlabs <- termlabs[-which_dynamics]
      for (i in seq_along(which_dynamics)) {
        fts_smooths[[i]] <- eval(parse(text = tf[which_dynamics[i]]))
      }
    }

    # Replace dynamic terms with the correct specification
    for (i in seq_along(which_dynamics)) {
      fts_forms <- dyn_to_spline(
        fts_smooths[[i]],
        term_id = i,
        data = data,
        formula = formula,
        time_var = time_var,
        gam_init = gam_init[[i]],
        newdata = newdata
      )
      gam_init[[i]] <- fts_forms$gam_init

      # Update formula
      termlabs <- c(termlabs, unlist(fts_forms$newform))

      # Update data
      if (inherits(data, "list")) {
        data <- append(data, as.list(fts_forms$X))
      } else {
        data <- dplyr::bind_cols(data, fts_forms$X)
      }
    }

    # Return the updated formula for passing to mgcv
    updated_formula <- reformulate(termlabs, rlang::f_lhs(formula))
    attr(updated_formula, ".Environment") <- attr(formula, ".Environment")
  } else {
    updated_formula <- formula
  }

  return(
    list(
      formula = updated_formula,
      data = data,
      fts_smooths = fts_smooths,
      gam_init = gam_init
    )
  )
}

#' Function to create the individual spline effects for time-varying
#' coefficients of basis functions
#' @noRd
dyn_to_spline <- function(
    term,
    term_id,
    data,
    formula,
    time_var = "time",
    gam_init = NULL,
    newdata = newdata) {
  # Extract key basis information
  label <- term$label
  time_k <- term$time_k
  time_bs <- term$time_bs

  # Initialise a gam object so the basis functions can be evaluated
  # and extracted; just use some Gaussian outcome here as all we need
  # is the predictor design matrix, with no column of 1s for intercepts
  if(is.null(newdata)) {
    data$my_fake_y <- rnorm(length(data[[1]]))
    gam_init <- ffc_gam_setup(
      formula(paste("my_fake_y ~ 0 + ", term$call)),
      dat = data
    )

    # Design matrix will contain the evaluated basis functions
    X <- predict(
      gam_init,
      type = "lpmatrix"
    )
  } else {
    X <- predict(
      gam_init,
      newdata = newdata,
      type = "lpmatrix"
    )
  }

  # Get indices of null space functions
  null_idx <- apply(X, 2, sd) == 0
  colnames(X) <- paste0(
    term$label,
    "_bs_",
    1:NCOL(X)
  )

  # Begin building the updated formula; start with any parametric
  # terms for null basis functions
  newform <- vector(mode = "list", length = NCOL(X))
  if (any(null_idx)) {
    for (i in which(null_idx)) {
      newform[[i]] <- colnames(X)[i]
    }
  }

  # Now build the smooth functions for time-varying coefficients
  for (i in which(!null_idx)) {
    newform[[i]] <- paste0(
      "s(",
      time_var,
      ", by = ",
      colnames(X)[i],
      ", ",
      "bs = '",
      time_bs,
      "', k = ",
      time_k,
      ", id = ",
      term_id,
      ")"
    )
  }

  # Return the formula list and the design matrix
  return(
    list(
      X = X,
      newform = newform,
      gam_init = gam_init
    )
  )
}
