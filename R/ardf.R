#' Fit an autoregressive dynamic factor model
#'
#' Fit an autoregressive dynamic factor model using Stan
#'
#' @importFrom fabletools new_model_class new_model_definition new_specials
#' @param formula Model specification (see "Specials" section)
#' @param family A family object specifying the outcome distribution to use in fitting.
#' Currently only `gaussian()` and `scat()` (i.e. Student-T) are supported
#' @param h `integer` specifying the forecast horizon
#' @param chains `integer` specifying the number of chains to be run
#' @param cores `integer` specifying the number of parallel cores to use
#' @param iter `integer` specifying the total number of iterations to run per chain
#' (including warmup)
#' @param warmup `integer` specifying the number of initial iterations to
#' use as burnin
#' @param adapt_delta the thin of the jumps in a HMC method
#' @param max_treedepth maximum tree depth per iteration
#' @param ... other arguments to pass to `rstan::sampling()`
#'
#' @author Nicholas J Clark
#'
#' @return A model specification
#' @examples
#' \donttest{
#' # Fit a functional forecasting model, then use ARDF for forecasting
#' library(dplyr)
#' 
#' # Split growth data into training and test sets
#' train_data <- growth_data %>% filter(age_yr <= 16)
#' test_data <- growth_data %>% filter(age_yr > 16)
#' 
#' # Step 1: Fit ffc_gam model with time-varying coefficients
#' mod <- ffc_gam(
#'   height_cm ~ fts(age_yr, by = id, time_k = 5),
#'   data = train_data,
#'   time = "age_yr",
#'   family = gaussian()
#' )
#' 
#' # Model fitted successfully
#' print("ffc_gam model completed")
#' 
#' # Step 2: Use ARDF for forecasting functional coefficients
#' # fc <- forecast(mod, newdata = test_data, model = "ARDF", 
#' #                chains = 1, iter = 500)
#' }
#' @export
ARDF = function(formula,
                family = gaussian(),
                h = 1,
                chains = 4,
                cores = 1,
                iter = 1000,
                warmup = floor(iter / 2),
                adapt_delta = 0.80,
                max_treedepth = 10,
                ...){

  ardf_model <- new_model_class(
    "ARDF",
    train = train_ardf,
    specials = specials_ardf,
    origin = NULL,
    check = all_tsbl_checks
  )
  new_model_definition(
    ardf_model,
    !!enquo(formula),
    family = family,
    h = h,
    chains = chains,
    cores = cores,
    iter = iter,
    warmup = warmup,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    ...
  )
}

#' @noRd
#' @return a tsibble containing forecast draws
#'
train_ardf = function(
    .data,
    specials,
    family = gaussian(),
    h = 1,
    chains = 4,
    cores = 1,
    iter = 1000,
    warmup = floor(iter / 2),
    adapt_delta = 0.80,
    max_treedepth = 10,
    ...){

  # Validate arguments
  checkmate::assert_count(chains, positive = TRUE)
  checkmate::assert_count(cores, positive = TRUE)
  checkmate::assert_count(iter, positive = TRUE)
  checkmate::assert_count(warmup, positive = TRUE)

  # Extract arguments from specials
  if(length(specials$K) > 1 || length(specials$p) > 1){
    rlang::warn("Only one special for `K()` and `p()` is allowed, defaulting to the first usage")
  }
  K <- specials$K[[1]]
  p <- specials$p[[1]]

  # Create the Stan data list
  model_data <- prep_tbl_ts_stan(
    .data = .data,
    h = h,
    K = K,
    p = p,
    family = family,
    model = 'ardf'
  )

  # Fit the model
  stanfit <- suppressWarnings(rstan::sampling(
    stanmodels$ardf,
    data = model_data,
    chains = chains,
    cores = cores,
    iter = iter,
    warmup = warmup,
    control = list(adapt_delta = adapt_delta,
                   max_treedepth = max_treedepth),
    ...
  ))

  # Extract forecasts as a tsibble
  out <- extract_stan_fc(
    stanfit,
    mod_name = 'ARDF',
    .data = .data,
    model_data = model_data,
    h = h
  )

  return(out)
}


#' @noRd
specials_ardf <- new_specials(
  K = function(K = 2) {
    checkmate::assert_count(K, positive = TRUE)
    as.list(environment())
  },
  p = function(p = 1) {
    checkmate::assert_count(p, positive = TRUE)
    as.list(environment())
  },
  .required_specials = c("K", "p")
)
