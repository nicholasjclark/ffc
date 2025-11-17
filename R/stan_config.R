#' Get Stan sampling defaults
#' @param type Character string specifying configuration type. Must be one of "sampling", "forecast", or "basis"
#' @return Named list of default parameter values
#' @examples
#' # Get sampling defaults
#' get_stan_defaults("sampling")
#' @noRd
get_stan_defaults <- function(type = "sampling") {
  checkmate::assert_choice(type, choices = c("sampling", "forecast", "basis"))
  
  switch(type,
    "sampling" = list(
      chains = 4,
      cores = 1,
      iter = 500,
      adapt_delta = 0.75,
      max_treedepth = 9
    ),
    "forecast" = list(
      h = 1,
      times = 200
    ),
    "basis" = list(
      time_k = 10
    )
  )
}

#' Get specific Stan parameter with fallback
#' @param param Character string specifying parameter name
#' @param type Character string specifying configuration type
#' @param default_value Fallback value if parameter not found
#' @return Parameter value or default_value if not found
#' @noRd
get_stan_param <- function(param, type = "sampling", default_value = NULL) {
  checkmate::assert_string(param)
  checkmate::assert_choice(type, choices = c("sampling", "forecast", "basis"))
  
  defaults <- get_stan_defaults(type)
  if (param %in% names(defaults)) {
    return(defaults[[param]])
  }
  default_value
}