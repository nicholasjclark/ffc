#'Extract model.frame from a fitted \code{ffc_gam} object
#'
#'
#'@inheritParams stats::model.frame
#'@param ... Ignored
#'@method model.frame ffc_gam
#'@author Nicholas J Clark
#'@return A \code{matrix} containing the fitted model frame
#'@export
model.frame.ffc_gam <- function(
    formula,
    ...) {

  # Remove class to make it easier to extract the original model.frame
  gam_object <- formula
  class(gam_object) <- class(gam_object)[-1]
  fts_smooths <- gam_object$fts_smooths

  if(is.null(fts_smooths)) {
    return(stats::model.frame(gam_object))

  } else {
    # If any fts() smooths included, extract their unique variable
    # names to ensure they are added appropriately to model.frame
    m_frame <- stats::model.frame(gam_object)
    all_terms <- unique(
      unlist(purrr::map(fts_smooths, 'term'))
    )

    all_terms <- unique(
      c(all_terms,
        unlist(purrr::map(fts_smooths, 'by')))
    )

    # Add any terms that aren't already in the model slot
    vars_to_add <- setdiff(
      all_terms,
      c('NA', colnames(m_frame))
    )

    if (length(vars_to_add)) {
      orig_names <- colnames(m_frame)
      for (i in 1:length(vars_to_add)) {
        gam_object$model <- cbind(
          m_frame,
          data[[vars_to_add[i]]]
        )
      }
      colnames(m_frame) <- c(orig_names, vars_to_add)
    }
    return(m_frame)
  }
}
