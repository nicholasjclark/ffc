#' Extract time-varying basis coefficients
#'
#' Extract time-varying basis coefficients a fitted \pkg{ffc_gam} model so they can
#' be visualized and / or forecasted ahead
#'
#' @name fts_coefs.ffc_gam
#' @importFrom mgcv predict.gam predict.bam
#' @param object \code{list} object of class \code{ffc_gam}. See [ffc_gam()]
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
#' @export
fts_coefs.ffc_gam <- function(
    object,
    ...) {
  if (is.null(object$fts_smooths)) {
    message("No functional smooths using fts() were included in this model")
    return(NULL)
  } else {
    # Time variable
    time_var <- object$time_var
    unique_times <- sort(unique(object$model[, time_var]))

    # Smooth names
    sm_names <- unlist(
      purrr::map(object$smooth, "label")
    )

    # fts() smooths
    fts_idx <- grep(":fts_",
      sm_names,
      fixed = TRUE
    )

    # Extract estimated coefficients
    betas <- mgcv::rmvn(
      n = 1000,
      mu = coef(object),
      V = vcov(object)
    )

    # Loop across fts smooths and calculate time-varying coefficients
    fts_preds <- do.call(
      rbind,
      lapply(fts_idx, function(sm) {
        # Create prediction grid
        by_var <- object$smooth[[sm]]$by
        pred_dat <- data.frame(
          by_var = 1,
          time_var = unique_times
        )
        colnames(pred_dat) <- c(by_var, time_var)

        # Create linear predictor matrix from the grid
        lp <- mgcv::PredictMat(
          object$smooth[[sm]],
          data = pred_dat,
        )

        # Calculate mean and SE of predictions
        beta_idx <- object$smooth[[sm]]$first.para:
        object$smooth[[sm]]$last.para
        preds <- matrix(
          NA,
          nrow = 1000,
          ncol = NROW(pred_dat)
        )
        for (i in 1:1000) {
          preds[i, ] <- as.vector(lp %*% betas[i, beta_idx])
        }

        # Return in tidy format
        dat <- data.frame(
          .basis = by_var,
          .time = unique_times,
          .estimate = apply(preds, 2, mean),
          .se = apply(preds, 2, function(x) sd(x) / sqrt(length(x)))
        )
        dat$time_var <- unique_times
        colnames(dat) <- c(
          ".basis",
          ".time",
          ".estimate",
          ".se",
          time_var
        )
        dat
      })
    )
    class(fts_preds) <- c("fts_ts", "tbl_df", "tbl", "data.frame")
    return(fts_preds)
  }
}

#' Print an fts_ts tibble
#'
#' @param x An object of class `fts_ts`
#' @param ... Ignored
#' @export
print.fts_ts <- function(x, ...) {
  NextMethod()
}

#' @export
format.fts_ts <- function(
    x, ...,
    n = NULL,
    width = NULL,
    n_extra = NULL) {
  NextMethod()
}

#' @noRd
fts_ts_2_tsbl <- function(x) {
  insight::check_if_installed("tsibble")
  tsibble::as_tsibble(
    functional_ts,
    key = .basis,
    index = .time
  )
}
