#' Fit a GP dymamic factor model
#'
#' Fit a squared exponential GP dymamic factor model using Stan
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
#' @export
GPDF = function(formula,
                family = gaussian(),
                h = 1,
                chains = 4,
                cores = 1,
                iter = 1000,
                warmup = floor(iter / 2),
                adapt_delta = 0.80,
                max_treedepth = 10,
                ...){

  gpdf_model <- new_model_class(
    "GPDF",
    train = train_gpdf,
    specials = specials_gpdf,
    origin = NULL,
    check = all_tsbl_checks
  )
  new_model_definition(
    gpdf_model,
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
train_gpdf = function(
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
  validate_pos_integer(chains)
  validate_pos_integer(cores)
  validate_pos_integer(iter)
  validate_pos_integer(warmup)

  # Extract arguments from specials
  if(length(specials$K) > 1){
    rlang::warn("Only one special for `K()` is allowed, defaulting to the first usage")
  }
  K <- specials$K[[1]]

  # Create the Stan data list
  model_data <- prep_tbl_ts_stan(
    .data = .data,
    h = h,
    K = K,
    p = 1,
    family = family,
    model = 'gpdf'
  )

  # Fit the model
  stanfit <- suppressWarnings(rstan::sampling(
    stanmodels$gpdf,
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
    mod_name = 'GPDF',
    .data = .data,
    model_data = model_data,
    h = h
  )

  return(out)
}


#' @noRd
specials_gpdf <- new_specials(
  K = function(K = 2) {
    validate_pos_integer(K)
    as.list(environment())
  },
  .required_specials = c("K")
)
