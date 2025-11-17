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

  .expanded <- .data |>
    tidyr::expand(tidyr::nesting(.basis, .realisation),
                  .time = min_time:max_time)
  .data |>
    as.data.frame() |>
    dplyr::right_join(.expanded,
                      by = dplyr::join_by(.basis, .time, .realisation)) |>
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
    out <- out |>
      dplyr::mutate(!!attr(.data, 'time_var') := .time)
  }
  return(out)
}

#' Checking functions written by Michell O'hara Wild and the fable team
#' @noRd
check_gaps <- function(x) {
  if (any(tsibble::has_gaps(x)[[".gaps"]])) {
    abort(sprintf("%s contains implicit gaps in time. You should check your data and convert implicit gaps into explicit missing values using `tsibble::fill_gaps()` if required.", deparse(substitute(x))))
  }
}

#' @noRd
check_regular <- function(x) {
  if (!tsibble::is_regular(x)) {
    abort(sprintf("%s is an irregular time series, which this model does not support. You should consider if your data can be made regular, and use `tsibble::update_tsibble(%s, regular = TRUE)` if appropriate.", deparse(substitute(x)), deparse(substitute(x))))
  }
}

#' @noRd
check_ordered <- function(x) {
  if (!tsibble::is_ordered(x)) {
    abort(sprintf(
      "%s is an unordered time series. To use this model, you first must sort the data in time order using `dplyr::arrange(%s, %s)`",
      deparse(substitute(x)), paste(c(deparse(substitute(x)), key_vars(x)), collapse = ", "), index_var(x)
    ))
  }
}

#' @noRd
all_tsbl_checks <- function(.data) {
  check_gaps(.data)
  check_regular(.data)
  check_ordered(.data)
  if (NROW(.data) == 0) {
    abort("There is no data to model. Please provide a dataset with at least one observation.")
  }
}

#' Extract response variable names from formula
#'
#' Handles both simple responses (y ~ ...) and cbind responses 
#' (cbind(successes, failures) ~ ...) robustly
#'
#' @param formula A formula object
#' @param return_all Logical, if TRUE returns all response variable names for 
#'   cbind, if FALSE reconstructs cbind() call (default FALSE)
#' @return Character vector of response variable name(s)
#' @noRd
extract_response_vars <- function(formula, return_all = FALSE) {
  # Identify response variable using robust approach
  resp_terms <- as.character(rlang::f_lhs(formula))
  
  if (length(resp_terms) == 1L) {
    response <- resp_terms
  } else {
    if (any(grepl('cbind', resp_terms))) {
      resp_terms <- resp_terms[-grepl('cbind', resp_terms)]
      
      if (return_all) {
        response <- resp_terms
      } else {
        # Reconstruct cbind() call as string
        response <- paste0('cbind(', paste(resp_terms, collapse = ','), ')')
      }
    } else {
      response <- resp_terms[1]
    }
  }
  
  return(response)
}
