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
  # Input validation
  if (is.list(formula)) {
    checkmate::assert_list(formula, types = "formula", min.len = 1,
                          .var.name = "formula")
  } else {
    checkmate::assert_class(formula, "formula", .var.name = "formula")
  }
  checkmate::assert_string(time_var, min.chars = 1, .var.name = "time_var")
  checkmate::assert_data_frame(data, min.rows = 1, .var.name = "data")
  checkmate::assert_list(gam_init, null.ok = TRUE, .var.name = "gam_init")
  if (!is.null(newdata)) {
    checkmate::assert_data_frame(newdata, .var.name = "newdata")
  }
  # Knots validation deferred to downstream ffc_gam_setup() following mgcv patterns
  
  # Save original data before adding row identifiers
  orig_data <- data
  
  # Add row identifiers for tracking through pipeline
  data <- add_row_identifiers(data)
  if (inherits(data, "tbl_df")) {
    data <- as.data.frame(data)
  }
  
  # Handle list formulae by processing each element
  if (is.list(formula)) {
    processed_formulae <- vector("list", length = length(formula))
    all_fts_smooths <- list()
    # Initialize parameter-aware structure for distributional models
    combined_gam_init <- vector("list", length = length(formula))
    
    # Ensure gam_init is properly structured for multiple formulae
    if (length(gam_init) == 0) {
      gam_init <- vector("list", length = length(formula))
    } else {
      # Validate existing structure for forecasting
      if (!validate_gam_init_structure(gam_init, length(formula))) {
        stop(insight::format_error(
          paste0("Invalid gam_init structure for distributional model with ",
                 length(formula), " parameters. Expected list of lists of GAM objects.")
        ), call. = FALSE)
      }
    }
    
    # Process each formula element
    for (i in seq_along(formula)) {
      current_formula <- formula[[i]]
      current_gam_init <- if (length(gam_init) >= i) {
        normalize_gam_init_structure(gam_init[[i]])
      } else {
        list()
      }
      
      # Process this formula element inline
      # Check if formula has an intercept
      keep_intercept <- attr(terms(current_formula), "intercept") == 1

      # Extract term labels
      termlabs <- attr(terms.formula(current_formula, keep.order = TRUE), 
                      "term.labels")

      # Check for offsets as well
      off_names <- grep(
        "offset",
        rownames(attr(terms.formula(current_formula), "factors")),
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
        for (j in seq_along(which_dynamics)) {
          fts_smooths[[j]] <- eval(parse(text = termlabs[which_dynamics[j]]))
        }

        # Replace dynamic terms with the correct specification
        fts_labs <- c()
        for (j in seq_along(which_dynamics)) {
          fts_forms <- dyn_to_spline(
            fts_smooths[[j]],
            term_id = j,
            data = data,
            formula = current_formula,
            time_var = time_var,
            gam_init = if (length(current_gam_init) >= j) {
              current_gam_init[[j]]
            } else {
              NULL
            },
            newdata = newdata,
            knots = knots,
            parameter_id = i
          )
          current_gam_init[[j]] <- fts_forms$gam_init

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
          rlang::f_lhs(current_formula)
        )
        attr(updated_formula, ".Environment") <- attr(current_formula, 
                                                     ".Environment")
      } else {
        updated_formula <- current_formula
      }

      single_result <- list(
        formula = updated_formula,
        data = data,
        fts_smooths = fts_smooths,
        gam_init = current_gam_init
      )
      
      # Store results
      processed_formulae[[i]] <- single_result$formula
      data <- single_result$data
      
      # Collect fts smooths with parameter indices
      if (length(single_result$fts_smooths) > 0) {
        param_smooths <- lapply(single_result$fts_smooths, function(smooth) {
          smooth$parameter <- i
          smooth
        })
        all_fts_smooths <- c(all_fts_smooths, param_smooths)
      }
      
      # Collect gam_init objects - maintain parameter structure for distributional models
      if (length(single_result$gam_init) > 0) {
        # Ensure combined_gam_init has proper parameter slots
        if (length(combined_gam_init) < i) {
          length(combined_gam_init) <- i
        }
        # Store gam_init objects for this parameter
        combined_gam_init[[i]] <- single_result$gam_init
      }
    }
    
    return(
      list(
        formula = processed_formulae,
        data = data,
        orig_data = orig_data,
        fts_smooths = all_fts_smooths,
        gam_init = combined_gam_init
      )
    )
  }

  # Process single formula (fallback for non-list case)
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
        gam_init = if (length(gam_init) >= i) gam_init[[i]] else NULL,
        newdata = newdata,
        knots = knots,
        parameter_id = NULL
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

#' Convert parameter ID to semantic name
#' @noRd
get_parameter_prefix <- function(parameter_id) {
  if (is.null(parameter_id)) {
    return(NULL)
  }
  
  semantic_names <- c("location", "scale", "shape")
  if (parameter_id <= length(semantic_names)) {
    return(semantic_names[parameter_id])
  } else {
    return(paste0("param", parameter_id))
  }
}

#' Normalize gam_init structure to ensure consistent list format
#' @param obj Object that could be a gam object, list, or NULL
#' @return List containing gam objects or empty list
#' @noRd
normalize_gam_init_structure <- function(obj) {
  checkmate::assert(
    checkmate::check_class(obj, "gam"),
    checkmate::check_list(obj),
    checkmate::check_null(obj)
  )
  
  if (inherits(obj, "gam")) {
    list(obj)
  } else if (is.list(obj)) {
    obj
  } else {
    list()
  }
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
    knots = NULL,
    parameter_id = NULL) {
  # Extract key basis information
  label <- term$label
  time_k <- term$time_k
  time_bs <- term$time_bs
  time_m <- term$time_m
  mean_only <- term$mean_only
  share_penalty <- term$share_penalty
  
  # For distributional families, default share_penalty to FALSE
  # Reason: shared penalties may cause fitting issues with mgcv 
  # distributional families
  if (!is.null(parameter_id) && share_penalty) {
    if (!identical(Sys.getenv("TESTTHAT"), "true")) {
      rlang::warn(
        "Shared penalties may cause fitting issues with distributional families. Setting {.field share_penalty} = FALSE.",
        .frequency = "once", 
        .frequency_id = "distributional_share_penalty"
      )
    }
    share_penalty <- FALSE
  }

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
  base_names <- paste0(
    "fts_bs_",
    clean_sm_names(names(coef(gam_init)))
  )
  
  # Add parameter prefix for distributional models
  param_prefix <- get_parameter_prefix(parameter_id)
  if (!is.null(param_prefix)) {
    base_names <- paste0(param_prefix, "_", base_names)
  }
  
  colnames(X) <- base_names

  if (mean_only) {
    # Ensure a constant basis is included to capture
    # any mean-shifts in the functions over time
    if (!is.null(newdata) || !any(apply(X, 2, sd) == 0)) {
      orig_names <- colnames(X)
      X <- cbind(
        X,
        matrix(1, nrow = NROW(X), ncol = 1)
      )
      mean_name <- paste0(label, term_id, "_mean")
      param_prefix <- get_parameter_prefix(parameter_id)
      if (!is.null(param_prefix)) {
        mean_name <- paste0(param_prefix, "_", mean_name)
      }
      
      colnames(X) <- c(orig_names, mean_name)
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
