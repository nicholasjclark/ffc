#' Fit a GP dynamic factor model
#'
#' Fit a squared exponential GP dynamic factor model using Stan
#'
#' @importFrom fabletools new_model_class new_model_definition new_specials
#' @inheritParams ARDF
#'
#' @author Nicholas J Clark
#'
#' @return A model specification
#' @examples
#' # Fit a functional forecasting model, then use GPDF for forecasting
#' library(dplyr)
#'
#' # Split growth data into training and test sets
#' train_data <- growth_data |> filter(age_yr <= 13)
#' test_data <- growth_data |> filter(age_yr > 13)
#'
#' # Step 1: Fit ffc_gam model with time-varying coefficients
#' mod <- ffc_gam(
#'   height_cm ~
#'     s(id, bs = "re") +
#'     fts(age_yr, time_k = 5),
#'   data = train_data,
#'   time = "age_yr",
#'   family = gaussian()
#' )
#'
#' # Step 2: Use GPDF for forecasting functional coefficients
#' fc <- forecast(mod, newdata = test_data, model = "GPDF",
#'                chains = 1, iter = 300)
#' @export
GPDF = function(
  formula,
  family = gaussian(),
  h = get_stan_param("h", "forecast"),
  chains = get_stan_param("chains"),
  cores = get_stan_param("cores"),
  iter = get_stan_param("iter"),
  warmup = floor(iter / 2),
  adapt_delta = get_stan_param("adapt_delta"),
  max_treedepth = get_stan_param("max_treedepth"),
  silent = get_stan_param("silent"),
  ...
) {
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
    silent = silent,
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
  h = get_stan_param("h", "forecast"),
  chains = get_stan_param("chains"),
  cores = get_stan_param("cores"),
  iter = get_stan_param("iter"),
  warmup = floor(get_stan_param("iter") / 2),
  adapt_delta = get_stan_param("adapt_delta"),
  max_treedepth = get_stan_param("max_treedepth"),
  silent = get_stan_param("silent"),
  ...
) {
  # Extract arguments from specials
  if (length(specials$K) > 1) {
    rlang::warn(
      "Only one special for `K()` is allowed, defaulting to the first usage"
    )
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
  stanfit <- run_stan_sampling(
    "gpdf",
    model_data,
    chains,
    cores,
    iter,
    warmup,
    adapt_delta,
    max_treedepth,
    silent,
    ...
  )

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
    checkmate::assert_count(K, positive = TRUE)
    as.list(environment())
  },
  .required_specials = c("K")
)
