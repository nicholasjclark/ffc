#' Forecasting functional basis coefficients
#'
#' @param object An object of class `fts_ts` containing time-varying
#' basis function coefficients extracted from an `ffc_gam` object
#' @param model A `character` string representing a valid univariate model definition
#' from the \pkg{fable} package or one of the built-in Bayesian dynamic
#' factor models. Note that if a \pkg{fable} model is used,
#' the chosen method must have an associated
#' `generate()` method in order to simulate forecast realisations. Valid models
#' currently include: `'ARDF'`, `'GPDF'`, '`VARDF`,
#' `'ETS'`, `'ARIMA'`, `'AR'`, `'RW'`, `'NAIVE'`, and `'NNETAR'`
#' @param h A positive `integer` specifying the length of the forecast
#' horizon
#' @param times A positive `integer` specifying the number of forecast
#' realisation paths to simulate from the fitted forecast `model`
#' @param stationary If `TRUE`, the fitted time series models are
#' constrained to be stationary. Default is `FALSE`. This option
#' only works when `model == 'ARIMA'`
#' @param ... Other arguments to pass to the Stan dynamic factor models
#' (i.e. the (V)AR order `lag = ...` or the number of factors `K = ...`)
#' @details Basis function coefficient time series will be used as input
#' to the specified `model` to train a forecasting model that will then be
#' extrapolated `h` timesteps into the future. A total of `times` forecast realisations
#' will be returned for each basis coefficient
#' @return A `tsibble` object containing the forecast prediction (`.sim`) for each
#' replicate realisation (`.rep`) at each timestep in the forecast horizon `h`
#' @seealso [predict()], [fts()]
#' @author Nicholas J Clark
#' @export
forecast.fts_ts <- function(
    object,
    model = "ARIMA",
    h = 1,
    times = 25,
    stationary = FALSE,
    ...
) {
  # Ensure required fable package is installed
  insight::check_if_installed("fable")

  # Validate forecast horizon and simulation count
  validate_pos_integer(h)
  validate_pos_integer(times)

  # Stationarity check for ARIMA model
  if (stationary & !identical(model, "ARIMA")) {
    stop("stationary = TRUE only works with model = 'ARIMA'")
  }

  # Convert functional time series to tsibble format
  object_tsbl <- fts_ts_2_tsbl(object)

  # Use Stan dynamic factor model for forecasts if specified
  if (model %in% c("ARDF", "VARDF", "GPDF")) {
    return(
      stan_forecast(object_tsbl, model, h, ...)
    )
  }

  # If using fable, check for uncertainty estimates and simulate forecasts
  if (exists(".sd", object_tsbl)) {
    object_sds <- object_tsbl %>%
      as.data.frame() %>%
      dplyr::group_by(.basis, .realisation) %>%
      dplyr::summarise(.sd = tail(.sd, 1))
  } else {
    object_sds <- NULL
  }

  return(
    fable_forecast(
      object_tsbl = object_tsbl,
      model = model,
      h = h,
      times = times,
      stationary = stationary,
      object_sds = object_sds
    )
  )
}

#' Extract model-specific parameters from function arguments
#' Defaults: K = 2 (number of factors), lag = 1 (time lag)
#' @noRd
extract_model_params <- function(...) {
  dots <- list(...)
  K <- if ("K" %in% names(dots)) dots$K else 2
  lag <- if ("lag" %in% names(dots)) dots$lag else 1
  list(
    K = K,
    p = lag
  )
}

#' Fit Stan dynamic factor model based on model type
#' Supported models: ARDF, VARDF, GPDF
#' @noRd
stan_forecast <- function(object_tsbl, model, h, ...) {
  params <- extract_model_params(...)

  if (model == "ARDF") {
    train_ardf(
      .data = object_tsbl,
      specials = params,
      h = h,
      family = gaussian(),
      chains = 4,
      cores = 4,
      iter = 500,
      warmup = 450,
      adapt_delta = 0.75,
      max_treedepth = 9
    )
  } else if (model == "VARDF") {
    train_vardf(
      .data = object_tsbl,
      specials = params,
      h = h,
      family = gaussian(),
      chains = 4,
      cores = 4,
      iter = 500,
      warmup = 450,
      adapt_delta = 0.75,
      max_treedepth = 9
    )
  } else if (model == "GPDF") {
    train_gpdf(
      .data = object_tsbl,
      specials = list(K = params$K),
      h = h,
      family = gaussian(),
      chains = 4,
      cores = 4,
      iter = 500,
      warmup = 450,
      adapt_delta = 0.75,
      max_treedepth = 9
    )
  }
}
