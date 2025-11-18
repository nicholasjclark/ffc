#' Interpret the formula specified to ffc and replace any fts() terms
#' with the correct time-varying coefficient terms
#' @importFrom stats formula terms as.formula terms.formula
#' @noRd
interpret_ffc <- function(
    formula,
    data,
    time_var = "time",
    gam_init = list(),
    newdata = NULL,
    knots = NULL) {
  # Tibbles often get rearranged by mgcv in very strange ways; grouping?
  # best to convert to plain-Jane data.frames first
  
  # Save original data before adding row identifiers
  orig_data <- data
  
  # Add row identifiers for tracking through pipeline
  data <- add_row_identifiers(data)
  if (inherits(data, "tbl_df")) {
    data <- as.data.frame(data)
  }

  # Check if formula has an intercept
  keep_intercept <- attr(terms(formula), "intercept") == 1

  # Extract term labels
  termlabs <- attr(terms.formula(formula, keep.order = TRUE), "term.labels")

  # Check for offsets as well
  off_names <- grep(
    "offset",
    rownames(attr(terms.formula(formula), "factors")),
    value = TRUE
  )
  if (length(off_names) > 0L) {
    termlabs <- c(termlabs, off_names)
  }

  # Check if any terms use the fts() wrapper
  which_dynamics <- grep("fts(", termlabs, fixed = TRUE)

  # Update the formula to the correct time-varying coefficient implementation
  fts_smooths <- list()
  if (length(which_dynamics) != 0L) {
    fts_smooths <- vector(length = length(which_dynamics), mode = "list")

    # Evaluate the fts() terms and update the formula appropriately
    for (i in seq_along(which_dynamics)) {
      fts_smooths[[i]] <- eval(parse(text = termlabs[which_dynamics[i]]))
    }

    # Replace dynamic terms with the correct specification
    fts_labs <- c()
    for (i in seq_along(which_dynamics)) {
      fts_forms <- dyn_to_spline(
        fts_smooths[[i]],
        term_id = i,
        data = data,
        formula = formula,
        time_var = time_var,
        gam_init = gam_init[[i]],
        newdata = newdata,
        knots = knots
      )
      gam_init[[i]] <- fts_forms$gam_init

      # Update formula
      fts_labs <- c(fts_labs, unlist(fts_forms$newform))

      # Update data
      if (inherits(data, "list")) {
        data <- append(data, as.list(fts_forms$X))
      } else {
        data <- cbind(data, fts_forms$X)
      }
    }

    # Return the updated formula for passing to mgcv
    updated_formula <- reformulate(
      c(termlabs[-which_dynamics], fts_labs),
      rlang::f_lhs(formula)
    )
    attr(updated_formula, ".Environment") <- attr(formula, ".Environment")
  } else {
    updated_formula <- formula
  }

  return(
    list(
      formula = updated_formula,
      data = data,
      orig_data = orig_data,
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
    newdata = newdata,
    knots = NULL) {
  # Extract key basis information
  label <- term$label
  time_k <- term$time_k
  time_bs <- term$time_bs
  time_m <- term$time_m
  mean_only <- term$mean_only
  share_penalty <- term$share_penalty

  # Initialise a gam object so the basis functions can be evaluated
  # and extracted; just use some Gaussian outcome here as all we need
  # is the predictor design matrix, with no column of 1s for intercepts
  if (is.null(newdata)) {
    data$my_fake_y <- rnorm(length(data[[1]]))
    myform <- formula(paste("my_fake_y ~ 0 + ", term$call),
      env = attr(formula, ".Environment")
    )
    gam_init <- ffc_gam_setup(
      myform,
      dat = data,
      knots = knots
    )

    # Design matrix will contain the evaluated basis functions
    X <- mgcv::predict.gam(
      gam_init,
      type = "lpmatrix"
    )
  } else {
    X <- mgcv::predict.gam(
      gam_init,
      newdata = newdata,
      type = "lpmatrix"
    )
  }
  
  # Validate basis matrix dimensions match data structure for row tracking
  expected_rows <- if (is.null(newdata)) nrow(data) else nrow(newdata)
  actual_rows <- nrow(X)
  
  if (actual_rows != expected_rows) {
    stop(insight::format_error(
      paste0("Basis matrix dimensions mismatch: got ", actual_rows, 
             " rows but expected ", expected_rows, 
             " rows to match data structure")
    ), call. = FALSE)
  }
  
  # Validate row count consistency for proper row ID mapping
  # Only warn if we expect 1:1 correspondence in training data
  if (!is.null(newdata) && ".original_row_id" %in% names(newdata)) {
    unique_ids <- length(unique(newdata$.original_row_id))
    # Only validate if expecting direct correspondence
    if (actual_rows != unique_ids && actual_rows == nrow(newdata)) {
      warning(insight::format_warning(
        paste0("Basis matrix rows (", actual_rows, 
               ") may not correspond to unique row IDs (", unique_ids, ")")
      ))
    }
  } else if (is.null(newdata) && ".row_id" %in% names(data)) {
    unique_ids <- length(unique(data$.row_id))
    if (actual_rows != unique_ids) {
      warning(insight::format_warning(
        paste0("Basis matrix rows (", actual_rows, 
               ") may not correspond to unique row IDs (", unique_ids, ")")
      ))
    }
  }

  # Get indices of null space functions
  colnames(X) <- paste0(
    "fts_bs_",
    clean_sm_names(names(coef(gam_init)))
  )

  if (mean_only) {
    # Ensure a constant basis is included to capture
    # any mean-shifts in the functions over time
    if (!is.null(newdata) || !any(apply(X, 2, sd) == 0)) {
      orig_names <- colnames(X)
      X <- cbind(
        X,
        matrix(1, nrow = NROW(X), ncol = 1)
      )
      colnames(X) <- c(
        orig_names,
        paste0(
          label,
          term_id,
          "_mean"
        )
      )
    }

    X <- X[, grepl(paste0(
      label,
      term_id,
      "_mean"
    ), colnames(X)),
    drop = FALSE
    ]
  }

  # Begin building the updated formula
  newform <- vector(mode = "list", length = NCOL(X))
  for (i in 1:NCOL(X)) {
    if (share_penalty) {
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
        ", m = ",
        time_m,
        ", id = ",
        term_id,
        ")"
      )
    } else {
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
        ", m = ",
        time_m,
        ")"
      )
    }
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

#' Clean smooth names so no illegal characters are used in formulae
#' @noRd
clean_sm_names <- function(sm_names) {
  sm_names_clean <- gsub(" ", "_", sm_names, fixed = TRUE)
  sm_names_clean <- gsub("(", "_", sm_names_clean, fixed = TRUE)
  sm_names_clean <- gsub(")", "_", sm_names_clean, fixed = TRUE)
  sm_names_clean <- gsub(",", "by", sm_names_clean, fixed = TRUE)
  sm_names_clean <- gsub(":", "by", sm_names_clean, fixed = TRUE)
  sm_names_clean <- gsub(".", "_", sm_names_clean, fixed = TRUE)
  sm_names_clean <- gsub("]", "_", sm_names_clean, fixed = TRUE)
  sm_names_clean <- gsub("[", "_", sm_names_clean, fixed = TRUE)
  sm_names_clean <- gsub(";", "_", sm_names_clean, fixed = TRUE)
  sm_names_clean <- gsub(":", "_", sm_names_clean, fixed = TRUE)
  sm_names_clean <- gsub("'", "", sm_names_clean, fixed = TRUE)
  sm_names_clean <- gsub("\"", "", sm_names_clean, fixed = TRUE)
  sm_names_clean <- gsub("%", "percent", sm_names_clean, fixed = TRUE)
  sm_names_clean <- gsub("[.]+", "_", sm_names_clean, fixed = TRUE)
  sm_names_clean <- gsub("'", "", sm_names_clean, fixed = TRUE)
  sm_names_clean
}
