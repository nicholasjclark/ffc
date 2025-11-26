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
      max_treedepth = 9,
      silent = TRUE
    ),
    "forecast" = list(
      h = 1,
      n_samples = 200
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

#' Execute Stan model sampling with standardized parameters
#' @param model_name Character string: model name ("ardf", "vardf", "gpdf")
#' @param model_data List: Stan data prepared by prep_tbl_ts_stan()
#' @param chains Integer: number of MCMC chains
#' @param cores Integer: number of CPU cores
#' @param iter Integer: total iterations per chain
#' @param warmup Integer: warmup iterations
#' @param adapt_delta Numeric: target acceptance rate
#' @param max_treedepth Integer: maximum tree depth
#' @param silent Logical: suppress progress output and viewer
#' @param ... Additional arguments passed to rstan::sampling()
#' @return rstan fit object
#' @noRd
run_stan_sampling <- function(model_name, model_data, chains, cores, iter, 
                              warmup, adapt_delta, max_treedepth, silent, ...) {
  # Comprehensive parameter validation
  checkmate::assert_choice(model_name, choices = c("ardf", "vardf", "gpdf"))
  checkmate::assert_list(model_data, min.len = 1)
  checkmate::assert_count(chains, positive = TRUE)
  checkmate::assert_count(cores, positive = TRUE)
  checkmate::assert_count(iter, positive = TRUE)
  checkmate::assert_count(warmup, positive = TRUE)
  checkmate::assert_number(adapt_delta, lower = 0, upper = 1)
  checkmate::assert_count(max_treedepth, positive = TRUE)
  checkmate::assert_logical(silent, len = 1)
  
  # Temporarily override rstan options and environment when silent
  if (silent) {
    old_messages <- getOption("rstan.show_messages", TRUE)
    options(rstan.show_messages = FALSE)
    on.exit(options(rstan.show_messages = old_messages), add = TRUE)
    
    # Temporarily disable RStudio detection by modifying environment
    old_rstudio <- Sys.getenv("RSTUDIO", unset = NA)
    if (!is.na(old_rstudio)) {
      Sys.unsetenv("RSTUDIO")
      on.exit(Sys.setenv(RSTUDIO = old_rstudio), add = TRUE)
    }
  }
  
  # Build sampling arguments
  sampling_args <- list(
    object = stanmodels[[model_name]],
    data = model_data,
    chains = chains,
    cores = cores,
    iter = iter,
    warmup = warmup,
    verbose = !silent,
    refresh = if (silent) 0 else max(1, iter %/% 10),
    show_messages = !silent,
    control = list(adapt_delta = adapt_delta,
                   max_treedepth = max_treedepth)
  )
  
  # Only set open_progress when we want to disable it
  # Let rstan handle progress automatically when silent = FALSE
  if (silent) {
    sampling_args$open_progress <- FALSE
  }
  
  # Add any additional arguments
  sampling_args <- c(sampling_args, list(...))
  
  suppressWarnings(do.call(rstan::sampling, sampling_args))
}