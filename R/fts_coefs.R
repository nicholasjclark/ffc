#' Extract parameter name from basis variable name
#'
#' @param basis_var Character string. Name of the basis variable
#' @param family mgcv family object to get semantic parameter names
#' @return Character string with parameter name (mu, sigma, etc.)
#' @noRd
extract_parameter_from_basis <- function(basis_var, family) {
  checkmate::assert_string(basis_var, min.chars = 1)

  # Standard parameter names: location, scale, shape
  standard_param_names <- c("location", "scale", "shape")

  # Check if this has a semantic parameter prefix (location_, scale_, shape_)
  for (param_name in standard_param_names) {
    if (grepl(paste0("^", param_name, "_"), basis_var)) {
      return(param_name)
    }
  }

  # Check if this has a numeric parameter prefix (paramN_)
  if (grepl("^param[0-9]+_", basis_var)) {
    # Extract parameter number
    param_num <- as.numeric(gsub("^param([0-9]+)_.*", "\\1", basis_var))

    # Use standard names based on parameter position
    if (param_num <= length(standard_param_names)) {
      return(standard_param_names[param_num])
    }

    # Fallback for parameters beyond shape
    return(paste0("param", param_num))
  }

  # Single-parameter model - always location
  return("location")
}

#' Extract time-varying basis coefficients
#'
#' Extract time-varying basis coefficients a fitted \pkg{ffc_gam} model so they can
#' be visualized and / or forecasted ahead
#'
#' @name fts_coefs.ffc_gam
#' @importFrom mgcv predict.gam predict.bam
#' @importFrom rlang := !! !!!
#' @param object \code{list} object of class \code{ffc_gam}. See [ffc_gam()]
#' @param summary `Logical`. Should summary statistics of the coefficient
#' time series be returned instead of realized curves? Default is `TRUE`. If
#' `FALSE`, replicate realisations of each basis coefficient time series will
#' be returned
#' @param n_samples A positive `integer` specifying the number of time series
#' realisation paths to simulate from the fitted model. Ignored if
#' `summary = FALSE`
#' @param ... Ignored
#' @details This function creates a `tidy` time series of basis function coefficients
#' from all `fts()` terms that were supplied to the original model
#' @return A `fts_ts` object containing the point estimates and their standard errors
#' for basis function coefficients
#' @seealso [ffc_gam()], [fts()]
#' @author Nicholas J Clark
#' @export
fts_coefs <- function(object, ...) {
  UseMethod("fts_coefs", object)
}

#' @rdname fts_coefs.ffc_gam
#' @method fts_coefs ffc_gam
#' @author Nicholas J Clark
#' @export
fts_coefs.ffc_gam <- function(
  object,
  summary = TRUE,
  n_samples = 25,
  ...
) {
  if (is.null(object$fts_smooths)) {
    message("No functional smooths using fts() were included in this model")
    return(NULL)
  } else {
    # Sanity check
    checkmate::assert_count(n_samples, positive = TRUE)

    # Time variable
    time_var <- object$time_var
    checkmate::assert_string(time_var, min.chars = 1)

    # Get all unique times within training window, assuming
    # equal time steps
    unique_times <- seq(
      min(object$model[[time_var]]),
      max(object$model[[time_var]])
    )

    # Extract unique index values for data that were provided as
    # tsibbles
    if (!is.null(attr(object$model, "index"))) {
      # Get the name of the index column
      index_name <- attr(object$model, "index")

      # Extract unique index values
      unique_indices <- unique(object$model[[index_name]])
    }

    # Smooth names
    sm_names <- unlist(
      purrr::map(object$smooth, "label")
    )

    # fts() smooths
    fts_idx <- grep("fts_", sm_names, fixed = TRUE)

    # Extract estimated coefficients
    betas <- mgcv::rmvn(
      n = n_samples,
      mu = coef(object),
      V = vcov(object)
    )

    # Use map_dfr() to iterate over each smooth index in fts_idx and
    # bind results row-wise
    fts_preds <- purrr::map_dfr(fts_idx, function(sm) {
      # Extract the 'by' variable name for the current smooth
      by_var <- object$smooth[[sm]]$by

      # Create a prediction grid with the 'by' variable set to 1 and
      # time variable set to unique time points
      pred_dat <- tibble::tibble(
        !!by_var := 1,
        !!time_var := unique_times
      )

      # Generate the linear predictor matrix for the current smooth
      lp <- mgcv::PredictMat(
        object$smooth[[sm]],
        data = pred_dat
      )

      # Identify the indices of the coefficients relevant to this smooth
      beta_idx <- object$smooth[[sm]]$first.para:object$smooth[[sm]]$last.para

      # Preallocate a matrix to store predictions for each simulation
      preds <- matrix(NA, nrow = n_samples, ncol = nrow(pred_dat))

      # Fill the prediction matrix by multiplying the predictor matrix with
      # the sampled coefficients
      for (i in 1:n_samples) {
        preds[i, ] <- as.vector(lp %*% betas[i, beta_idx])
      }

      # Extract parameter name from basis variable name
      parameter_name <- extract_parameter_from_basis(by_var, object$family)

      # If summary output is requested, compute mean and standard error
      # across simulations
      if (summary) {
        dat <- tibble::tibble(
          .basis = by_var,
          .parameter = parameter_name,
          .time = unique_times,
          .estimate = apply(preds, 2, mean),
          .se = apply(preds, 2, function(x) sd(x) / sqrt(length(x))),
          !!time_var := unique_times
        )
      } else {
        # Otherwise, return all simulation replicates in long format
        # using map_dfr()
        dat <- purrr::map_dfr(1:n_samples, function(x) {
          repdat <- tibble::tibble(
            .basis = by_var,
            .parameter = parameter_name,
            .time = unique_times,
            .estimate = preds[x, ],
            .realisation = x,
            !!time_var := unique_times
          )

          # If the model has an index attribute, add it to the output
          if (!is.null(attr(object$model, "index"))) {
            repdat <- dplyr::mutate(
              repdat,
              !!index_name := unique_indices
            )
          }

          repdat
        })
      }

      # Return the tidy data for this smooth
      dat
    })

    # Structure the object and return
    fts_preds <- structure(
      fts_preds,
      class = c("fts_ts", "tbl_df", "tbl", "data.frame"),
      time_var = time_var,
      index = attr(object$model, "index"),
      index2 = attr(object$model, "index2"),
      interval = attr(object$model, "interval"),
      summarized = isTRUE(summary)
    )

    return(fts_preds)
  }
}

#' Print an fts_ts tibble
#'
#' @param x An object of class `fts_ts`
#' @param ... Ignored
#' @author Nicholas J Clark
#' @export
print.fts_ts <- function(x, ...) {
  NextMethod()
}

#' @author Nicholas J Clark
#' @export
format.fts_ts <- function(
  x,
  ...,
  n = NULL,
  width = NULL,
  n_extra = NULL
) {
  NextMethod()
}

#' @noRd
fts_ts_2_tsbl <- function(x) {
  insight::check_if_installed("tsibble")
  if (attr(x, "summarized")) {
    out <- tsibble::as_tsibble(
      x,
      key = .basis,
      index = .time
    )
  } else {
    out <- tsibble::as_tsibble(
      x,
      key = c(.basis, .realisation),
      index = .time
    )
  }

  if (!is.null(attr(x, "interval"))) {
    attr(out, "interval") <- attr(x, "interval")
  }

  if (!is.null(attr(x, "index"))) {
    attr(out, "index") <- attr(x, "index")
  }

  if (!is.null(attr(x, "index2"))) {
    attr(out, "index2") <- attr(x, "index2")
  }

  return(out)
}
