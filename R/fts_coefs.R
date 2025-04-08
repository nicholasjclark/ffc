#' Extract time-varying basis coefficients
#'
#' Extract time-varying basis coefficients a fitted \pkg{ffc_gam} model so they can
#' be visualized and / or forecasted ahead
#'
#' @name fts_coefs.ffc_gam
#' @importFrom mgcv predict.gam predict.bam
#' @param object \code{list} object of class \code{ffc_gam}. See [ffc_gam()]
#' @param se.fit \code{logical} indicating whether standard errors of time-varying
#' basis coefficients should also be returned. Default is `TRUE`
#' @param ... Ignored
#' @details This function returns predictions from models fitted with [ffc_gam()].
#' Data passed to `newdata` will first be correctly augmented to include any basis functions
#' whose coefficients were estimated as time-varying, so the user need only supply data
#' that includes variables that were used in the original `data` that was supplied to
#' `ffc_gam()`
#' @return If `type == "lpmatrix"` then a `matrix` is returned which
#' will give a vector of linear predictor values (minus any offest)
#' at the supplied covariate values, when applied to the model coefficient vector.
#' Otherwise, if `se.fit == TRUE` then a `2 `item `list` is returned with
#' items (both `arrays`) `fit` and `se.fit` containing predictions and
#' associated standard error estimates, otherwise an `array` of predictions is
#' returned. The dimensions of the returned `arrays` depends on whether
#' type is `"terms"` or not: if it is then the `array` is `2` dimensional with each
#' term in the linear predictor separate, otherwise the `array` is 1 dimensional and
#' contains the linear predictor/predicted values (or corresponding s.e.s).
#' The linear predictor returned termwise will not include the offset or the intercept.
#' @seealso [ffc_gam()], [fts()]
#' @author Nicholas J Clark
#' @export
fts_coefs <- function(object, ...) {
  UseMethod("fts_coefs", object)
}

#' @rdname fts_coefs.ffc_gam
#' @method fts_coefs ffc_gam
#' @export
fts_coefs.ffc_gam = function(
    object,
    se.fit = TRUE,
    ...){

  if(is.null(object$fts_smooths)){
    message('No functional smooths using fts() were included in this model')
    return(NULL)
  } else {

    # Time variable
    time_var <- object$time_var
    unique_times <- sort(unique(object$model[, time_var]))

    # Smooth names
    sm_names <- unlist(
      purrr::map(object$smooth, 'label')
    )

    # fts() smooths
    fts_idx <- grep(':fts_',
                    sm_names,
                    fixed = TRUE)

    # Extract estimated coefficients
    betas <- mgcv::rmvn(
      n = 500,
      mu = coef(object),
      V = vcov(object)
    )

    # Loop across fts smooths and calculate time-varying coefficients
    fts_preds <- do.call(
      rbind,
      lapply(fts_idx, function(sm){

        # Create prediction grid
        by_var <- object$smooth[[sm]]$by
        pred_dat <- data.frame(by_var = 1,
                               time_var = unique_times)
        colnames(pred_dat) <- c(by_var, time_var)

        # Create linear predictor matrix from the grid
        lp <- mgcv::PredictMat(
          object$smooth[[sm]],
          data = pred_dat,
        )

        # Calculate mean and SE of predictions
        beta_idx <- object$smooth[[sm]]$first.para :
          object$smooth[[sm]]$last.para
        preds <- matrix(NA,
                        nrow = 500,
                        ncol = NROW(pred_dat))
        for (i in 1:500) {
          preds[i, ] <- as.vector(lp %*% betas[i, beta_idx])
        }

        # Return in tidy format
        data.frame(
          .basis = by_var,
          .time = unique_times,
          .estimate = apply(preds, 2, mean),
          .se = apply(preds, 2, function(x) sd(x) / length(x))
        )
      })
    )
    class(fts_preds) <- c('tbl_df', 'tbl', 'data.frame')
    return(fts_preds)
  }
}
