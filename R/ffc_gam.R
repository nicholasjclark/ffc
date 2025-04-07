#' Define functions with dynamic coefficients in \pkg{ffc} formulae
#'
#' Set up smooth terms with time-varying (dynamic) coefficients for use in
#' \pkg{ffc} models
#'
#' @importFrom mgcv gam bam
#' @inheritParams mgcv::bam
#' @param formula A GAM formula (see \code{\link[mgcv]{formula.gam}} and also
#' \code{\link[mgcv]{gam.models}}). This is exactly like the formula for a
#' GLM except that smooth terms, [fts()], \code{\link[mgcv]{s}()} and
#' \code{\link[mgcv]{te}()} can be added to the right hand side to specify
#' that the linear predictor depends on smooth functions of predictors
#' (or linear functionals of these).
#' @param engine `character` string specifying which \pkg{mgcv} interface to use
#' for fitting the model.
#' @param time `character` specifying which variable in `data` represents the
#' the time ordering of the observations
#' @param ... 	other arguments to pass to either \code{\link[mgcv]{gam}}
#' @rdname ffc_gam
#' @details This function will update the supplied `formula` to ensure any time-varying
#' functionals (supplied through [fts()] terms in the formula right hand side) are
#' appropriately incorporated into the model. It then passes the updated model and data
#' objects to the specified `engine` for model fitting
#' @return An object of class `ffc_gam`
#' @seealso [fts()]
#' @author Nicholas J Clark
#' @export
ffc_gam <- function(
    formula,
    family = gaussian(),
    data = list(),
    time,
    engine = c('gam', 'bam'),
    ...){

  engine <- match.arg(engine)
  orig_call <- match.call()
  orig_formula <- formula
  if(!exists(time, data)){
    stop(paste0(
      "the variable '",
      time,
      "' cannot be found in data"
    ),
    call. = FALSE)
  }

  # Update formula and data by checking for any fts() terms
  interpreted <- interpret_ffc(
    formula = formula,
    data = data,
    time_var = time
  )

  # Fit the model using the specified engine
  dots <- list(...)
  fit_args <- list(
    formula = interpreted$formula,
    family = family,
    data = interpreted$data,
    ...
  )
  out <- do.call(engine, fit_args)

  # Update the object and return
  out$call <- orig_call
  out$orig_formula <- orig_formula
  out$time_var <- time
  out$fts_smooths <- interpreted$fts_smooths
  out$gam_init <- interpreted$gam_init
  class(out) <- c('ffc_gam', class(out))
  return(out)
}
