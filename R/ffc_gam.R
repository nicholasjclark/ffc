#' Fit a functional time series model using dynamic functional coefficients
#'
#' Fit Generalized Additive Models that can include time-varying (dynamic)
#' functions
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
#' @return An object of class `ffc_gam`, which inherits from objects of class `gam` or
#' `bam`
#' @seealso [fts()], \code{\link[mgcv]{gam}}, \code{\link[mgcv]{bam}}
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
    na.action = 'na.fail',
    ...
  )
  out <- do.call(engine, fit_args)

  # Update the object and return
  out$call <- orig_call
  out$orig_formula <- orig_formula
  out$time_var <- time
  out$fts_smooths <- interpreted$fts_smooths
  out$gam_init <- interpreted$gam_init
  out <- update_mod_data(
    gam_object = out,
    fts_smooths = interpreted$fts_smooths,
    data = data
  )

  class(out) <- c('ffc_gam', class(out))
  return(out)
}

#' Ensure all terms are included in the stored model data
#' @noRd
update_mod_data <- function(
    gam_object,
    fts_smooths,
    data) {

  # All terms included in fts smooths
  all_terms <- unique(
    unlist(purrr::map(fts_smooths, 'term'))
  )

  all_terms <- unique(
    c(all_terms,
      unlist(purrr::map(fts_smooths, 'by')))
  )

  # Any offset terms
  termlabs <- attr(terms.formula(gam_object$formula, keep.order = TRUE), "term.labels")

  # Check for offsets as well
  off_names <- grep(
    'offset',
    rownames(attr(terms.formula(gam_object$formula), 'factors')),
    value = TRUE
  )
  if (length(off_names) > 0L) {
    all_terms <- c(all_terms, strip_offset(off_names))
  }

  # Add any terms that aren't already in the model slot
  vars_to_add <- setdiff(
    all_terms,
    c('NA', colnames(gam_object$model))
  )

  if (length(vars_to_add)) {
    orig_names <- colnames(gam_object$model)
    for (i in 1:length(vars_to_add)) {
      gam_object$model <- cbind(
        gam_object$model,
        data[[vars_to_add[i]]]
      )
    }
    colnames(gam_object$model) <- c(orig_names, vars_to_add)
  }
  return(gam_object)
}

#' Strip offset names down to the actual terms
#' @noRd
strip_offset <- function(x) {
  for (i in 1:length(x)) {
    if (substr(x[i], 1, 7) == "offset(")
      x[i] <- substr(x[i], 8, nchar(x[i]) - 1)

    if (substr(x[i], 1, 4) == "log(")
      x[i] <- substr(x[i], 5, nchar(x[i]) - 1)
  }
  x
}
