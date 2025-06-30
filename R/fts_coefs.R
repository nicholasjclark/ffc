#' Extract time-varying basis coefficients
#'
#' Extract time-varying basis coefficients a fitted \pkg{ffc_gam} model so they can
#' be visualized and / or forecasted ahead
#'
#' @name fts_coefs.ffc_gam
#' @importFrom mgcv predict.gam predict.bam
#' @importFrom rlang := !!
#' @param object \code{list} object of class \code{ffc_gam}. See [ffc_gam()]
#' @param summary `Logical`. Should summary statistics of the coefficient
#' time series be returned instead of realized curves? Default is `TRUE`. If
#' `FALSE`, replicate realisations of each basis coefficient time series will
#' be returned
#' @param times A positive `integer` specifying the number of time series
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
    times = 25,
    ...) {
  if (is.null(object$fts_smooths)) {
    message("No functional smooths using fts() were included in this model")
    return(NULL)
  } else {
    # Time variable
    time_var <- object$time_var

    # Get all unique times within training window, assuming
    # equal time steps
    unique_times <- seq(
      min(object$model[, time_var]),
      max(object$model[, time_var])
    )

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
      n = 500,
      mu = coef(object),
      V = vcov(object)
    )

    # Loop across fts smooths and calculate time-varying coefficients
    validate_pos_integer(times)
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
          nrow = 500,
          ncol = NROW(pred_dat)
        )
        for (i in 1:500) {
          preds[i, ] <- as.vector(lp %*% betas[i, beta_idx])
        }

        # Return in tidy format
        if (summary) {
          # dat <- data.frame(
          #   .basis = by_var,
          #   .time = unique_times,
          #   .estimate = apply(preds, 2, mean),
          #   .se = apply(preds, 2, function(x) sd(x) / sqrt(length(x)))
          # )
          # dat$time_var <- unique_times
          # colnames(dat) <- c(
          #   ".basis",
          #   ".time",
          #   ".estimate",
          #   ".se",
          #   time_var
          # )

          # Create a tibble summarizing the predictions across replications
          dat <- tibble::tibble(
            # Assign the basis variable (e.g., spline basis or grouping variable)
            .basis = by_var,

            # Assign the unique time points
            .time = unique_times,

            # Compute the mean estimate across replications for each time point (column-wise mean)
            .estimate = apply(preds, 2, mean),

            # Compute the standard error for each time point:
            # standard deviation divided by the square root of the number of replications
            .se = apply(preds, 2, function(x) sd(x) / sqrt(length(x))),

            # Dynamically name the time variable column using `:=`
            # and `!!` from `rlang`
            !!time_var := unique_times
          )

        } else {
          # dat <- do.call(
          #   rbind,
          #   lapply(1:times, function(x) {
          #     repdat <- data.frame(
          #       .basis = by_var,
          #       .time = unique_times,
          #       .estimate = preds[x, ],
          #       .rep = x
          #     )
          #     repdat$time_var <- unique_times
          #     colnames(repdat) <- c(
          #       ".basis",
          #       ".time",
          #       ".estimate",
          #       ".realisation",
          #       time_var
          #     )
          #     if (!is.null(attr(object$model, "index"))) {
          #       repdat[[attr(object$model, "index")]] <-
          #         unique(object$model[[attr(object$model, "index")]])
          #     }
          #
          #     repdat
          #   })
          # )

          # Use map_dfr() to iterate over 1:times and bind the resulting data
          # frames row-wise
          dat <- purrr::map_dfr(1:times, function(x) {

            # Create a tibble for the current replication
            # .basis, .time, and .estimate are populated from vectors
            # .realisation is the current iteration index
            # `!!time_var := unique_times` dynamically names the column using
            # the value of `time_var`
            repdat <- tibble::tibble(
              .basis = by_var,
              .time = unique_times,
              .estimate = preds[x, ],
              .realisation = x,
              !!time_var := unique_times
            )

            # If the model has an "index" attribute, add a column with
            # its unique values
            # This is useful for time series or panel data where an index
            # (e.g., subject ID) is needed
            if (!is.null(attr(object$model, "index"))) {

              # Get the name of the index column
              index_name <- attr(object$model, "index")

              # Add the index values
              repdat <- dplyr::mutate(
                repdat,
                !!index_name := unique(object$model[[index_name]])
              )
              #repdat[[index_name]] <- unique(object$model[[index_name]])
            }

            # Return the tibble for this iteration
            repdat
          })
        }
        dat
      })
    )
    class(fts_preds) <- c("fts_ts", "tbl_df", "tbl", "data.frame")
    attr(fts_preds, "time_var") <- time_var
    attr(fts_preds, "index") <- attr(object$model, "index")
    attr(fts_preds, "index2") <- attr(object$model, "index2")
    attr(fts_preds, "interval") <- attr(object$model, "interval")
    attr(fts_preds, "summarized") <- isTRUE(summary)
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
    x, ...,
    n = NULL,
    width = NULL,
    n_extra = NULL) {
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
