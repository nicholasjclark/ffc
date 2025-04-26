#' Predict from a fitted \pkg{ffc_gam} model
#'
#'
#' @importFrom mgcv predict.gam predict.bam
#' @importFrom stats predict
#' @inheritParams mgcv::predict.gam
#' @param ... ignored
#' @rdname predict.ffc_gam
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
#' @method predict ffc_gam
#' @author Nicholas J Clark
#' @export
predict.ffc_gam <- function(
    object,
    newdata,
    type = "link",
    se.fit = FALSE,
    ...) {
  type <- match.arg(
    type,
    choices = c(
      "link",
      "terms",
      "response",
      "lpmatrix",
      "iterms"
    )
  )

  # Update the supplied newdata in light of any fts() smooths
  if (missing(newdata)) {
    interpreted <- list()
    interpreted$data <- object$model
  } else {
    interpreted <- interpret_ffc(
      formula = object$orig_formula,
      data = newdata,
      newdata = newdata,
      gam_init = object$gam_init,
      time_var = object$time_var
    )
  }

  # Determine the prediction engine
  if (class(object)[2] == "gam") {
    pred_engine <- "predict.gam"
  } else {
    pred_engine <- "predict.bam"
  }

  # Compute predictions
  out <- do.call(
    pred_engine,
    list(
      object = object,
      newdata = interpreted$data,
      type = type,
      se.fit = se.fit
    )
  )

  return(out)
}
