#' Create distribution object from predictions
#' @param predictions Matrix of predictions
#' @return Distribution object
#' @noRd
create_distribution_object <- function(predictions) {
  distributional::dist_sample(
    lapply(seq_len(ncol(predictions)), function(i) predictions[, i])
  )
}

#' Compute summary statistics
#' @param predictions Matrix of predictions
#' @param robust Whether to use robust statistics
#' @return List with estimates and errors
#' @noRd
compute_summary_statistics <- function(predictions, robust) {
  if (robust) {
    estimates <- apply(predictions, 2, median, na.rm = TRUE)
    errors <- apply(abs(predictions - estimates), 2, median, na.rm = TRUE)
  } else {
    estimates <- apply(predictions, 2, mean, na.rm = TRUE)
    errors <- apply(predictions, 2, sd, na.rm = TRUE)
  }

  list(estimates = estimates, errors = errors)
}

#' Compute quantiles
#' @param predictions Matrix of predictions
#' @param probs Quantile probabilities
#' @return Matrix of quantiles
#' @noRd
compute_quantiles <- function(predictions, probs) {
  apply(predictions, 2, quantile, probs = probs, na.rm = TRUE)
}

#' Create summary data frame
#' @param stats Summary statistics
#' @param quantiles Quantile matrix
#' @param probs Quantile probabilities
#' @return Formatted data frame
#' @noRd
create_summary_dataframe <- function(stats, quantiles, probs) {
  out <- data.frame(
    .estimate = stats$estimates,
    .error = stats$errors,
    t(quantiles)
  )

  # Set column names
  colnames(out) <- c(
    ".estimate",
    ".error",
    paste0(".q", formatC(100 * probs, format = "f", digits = 1))
  )

  # Convert to tibble
  structure(out, class = c("tbl_df", "tbl", "data.frame"))
}

#' Format forecast output according to user preferences
#'
#' @param predictions Matrix of predictions
#' @param summary Whether to return summary statistics
#' @param robust Whether to use robust statistics
#' @param probs Quantile probabilities
#' @return Formatted forecast output
#' @noRd
format_forecast_output <- function(predictions,
                                   summary = TRUE,
                                   robust = FALSE,
                                   probs = c(0.025, 0.1, 0.9, 0.975)) {

  if (!summary) {
    return(create_distribution_object(predictions))
  }

  # Compute summary statistics
  stats <- compute_summary_statistics(predictions, robust)
  quantiles <- compute_quantiles(predictions, probs)

  # Create output data frame
  create_summary_dataframe(stats, quantiles, probs)
}

#' Get linear predictor matrix from model
#' @param object Model object
#' @param newdata New data
#' @return Linear predictor matrix
#' @noRd
get_linear_predictor_matrix <- function(object, newdata) {
  predict(object, newdata = newdata, type = "lpmatrix")
}

#' Preserve FTS attributes
#' @param functional_coefs New functional coefficients
#' @param intermed_coefs Original intermediate coefficients
#' @return Functional coefficients with preserved attributes
#' @noRd
preserve_fts_attributes <- function(functional_coefs, intermed_coefs) {
  # Add index information if needed
  if (!is.null(attr(intermed_coefs, "index"))) {
    if (!attr(intermed_coefs, "index") %in% names(functional_coefs)) {
      functional_coefs <- functional_coefs %>%
        dplyr::left_join(
          intermed_coefs %>%
            dplyr::select(.time, !!attr(intermed_coefs, "index")) %>%
            dplyr::distinct(),
          by = dplyr::join_by(.time)
        )
    }
  }

  # Set class and attributes
  structure(
    functional_coefs,
    class = c("fts_ts", "tbl_df", "tbl", "data.frame"),
    time_var = attr(intermed_coefs, "time_var"),
    index = attr(intermed_coefs, "index"),
    index2 = attr(intermed_coefs, "index2"),
    interval = attr(intermed_coefs, "interval"),
    summarized = attr(intermed_coefs, "summarized")
  )
}

#' Normalise quantiles into sampling weights
#' @noRd
norm_quantiles <- function(x) {
  xhat <- vector(length = length(x))
  for (i in seq_along(x)) {
    if (x[i] < 50) {
      xhat[i] <- (50 - x[i]) / 50
    } else {
      xhat[i] <- (x[i] - 50) / 50
    }
  }
  1.1 - xhat
}

#' @noRd
nlist = function(...) {
  m <- match.call()
  dots <- list(...)
  no_names <- is.null(names(dots))
  has_name <- if (no_names) FALSE else nzchar(names(dots))
  if (all(has_name)) return(dots)
  nms <- as.character(m)[-1]
  if (no_names) {
    names(dots) <- nms
  } else {
    names(dots)[!has_name] <- nms[!has_name]
  }
  dots
}

#' @noRd
`c<-` = function(x, value) {
  c(x, value)
}

#' Expand a tbl_ts to include out of sample rows with NAs
#' @noRd
expand_tbl_ts = function(.data, h){
  min_time <- min(.data$.time)
  max_time <- max(.data$.time) + h

  .expanded <- .data %>%
    tidyr::expand(tidyr::nesting(.basis, .realisation),
                  .time = min_time:max_time)
  .data %>%
    as.data.frame() %>%
    dplyr::right_join(.expanded,
                      by = dplyr::join_by(.basis, .time, .realisation)) %>%
    dplyr::ungroup()
}

#' Expand a tbl_ts to include out of sample index variables
#' @noRd
make_future_data <- function(.data, h = NULL){
  n <- fabletools::get_frequencies(h, .data, .auto = "smallest")
  if(length(n) > 1){
    warn("More than one forecast horizon specified, using the smallest.")
    n <- min(n)
  }
  if(is.null(h)) n <- n*2

  out <- tsibble::new_data(.data, round(n))
  if (tsibble::index_var(.data) == '.time') {
    out <- out %>%
      dplyr::mutate(!!attr(.data, 'time_var') := .time)
  }
  return(out)
}

#' Add model offset to linear predictors
#' @param linpreds Linear predictors matrix
#' @param linpred_matrix Original linear predictor matrix with offset attribute
#' @return Modified linear predictors matrix
#' @noRd
add_model_offset <- function(linpreds, linpred_matrix) {
  model_offset <- attr(linpred_matrix, "model.offset")

  if (is.null(model_offset)) {
    return(linpreds)
  }

  if (length(model_offset) == 1) {
    # Single offset value - add to all predictions
    return(linpreds + model_offset)
  } else if (length(model_offset) == ncol(linpreds)) {
    # Vector of offsets - add to each column
    return(sweep(linpreds, 2, model_offset, "+"))
  } else {
    warning("Model offset length (", length(model_offset),
            ") does not match number of predictions (", ncol(linpreds), ")")
    return(linpreds)
  }
}
