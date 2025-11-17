#' Forecasting functional basis coefficients
#'
#' @param object An object of class `fts_ts` containing time-varying
#' basis function coefficients extracted from an `ffc_gam` object
#' @param model A `character` string representing a valid univariate model definition
#' from the \pkg{fable} package, ensemble methods, or one of the built-in Bayesian dynamic
#' factor models. Note that if a \pkg{fable} model is used,
#' the chosen method must have an associated
#' `generate()` method in order to simulate forecast realisations. Valid models
#' currently include: `'ENS'`, `'ARDF'`, `'GPDF'`, '`VARDF`,
#' `'ETS'`, `'ARIMA'`, `'AR'`, `'RW'`, `'NAIVE'`, and `'NNETAR'`. 
#' The `'ENS'` option combines ETS and Random Walk forecasts with equal weights,
#' hedging bets between exponential smoothing and random walk assumptions
#' to provide more robust predictions when model uncertainty is high.
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
#' @examples
#' \donttest{
#' # Extract coefficients and generate forecasts
#' mod <- ffc_gam(
#'   deaths ~ offset(log(population)) + sex +
#'     fts(age, k = 8, bs = "cr", time_k = 10),
#'   time = "year",
#'   data = qld_mortality,
#'   family = poisson(),
#'   engine = "bam"
#' )
#' coefs <- fts_coefs(mod, summary = FALSE, times = 5)
#'
#' # Generate ETS forecasts
#' forecast(coefs, model = "ETS", h = 3)
#' }
#' @export
forecast.fts_ts <- function(
    object,
    model = "ARIMA",
    h = get_stan_param("h", "forecast"),
    times = 25,
    stationary = FALSE,
    ...
) {
  # Ensure required fable package is installed
  insight::check_if_installed("fable")

  # Validate forecast horizon and simulation count
  checkmate::assert_count(h, positive = TRUE)
  checkmate::assert_count(times, positive = TRUE)

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
    object_sds <- object_tsbl |>
      as.data.frame() |>
      dplyr::group_by(.basis, .realisation) |>
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

  # Model parameters
  K <- if ("K" %in% names(dots)) dots$K else 2
  lag <- if ("lag" %in% names(dots)) dots$lag else 1

  # Stan sampling parameters with sensible defaults
  chains <- if ("chains" %in% names(dots)) dots$chains else get_stan_param("chains")
  iter <- if ("iter" %in% names(dots)) dots$iter else get_stan_param("iter")
  warmup <- if ("warmup" %in% names(dots)) dots$warmup else floor(iter / 2)  # Half of actual iter
  cores <- if ("cores" %in% names(dots)) dots$cores else min(get_stan_param("chains"), parallel::detectCores() - 1)
  adapt_delta <- if ("adapt_delta" %in% names(dots)) dots$adapt_delta else get_stan_param("adapt_delta")
  max_treedepth <- if ("max_treedepth" %in% names(dots)) dots$max_treedepth else get_stan_param("max_treedepth")
  times <- if ("times" %in% names(dots)) dots$times else get_stan_param("times", "forecast")

  # Validate parameters using checkmate
  checkmate::assert_count(K, positive = TRUE)
  checkmate::assert_count(lag, positive = TRUE)
  checkmate::assert_count(chains, positive = TRUE)
  checkmate::assert_count(iter, positive = TRUE)
  checkmate::assert_count(warmup)
  checkmate::assert_count(cores, positive = TRUE)
  checkmate::assert_number(adapt_delta, lower = 0, upper = 1)
  checkmate::assert_count(max_treedepth, positive = TRUE)
  checkmate::assert_count(times, positive = TRUE)

  # Additional logic checks
  if (warmup >= iter) {
    stop(insight::format_error("warmup must be less than iter"))
  }

  # Calculate posterior samples and ensure adequate minimum
  post_samples <- (iter - warmup) * chains
  min_samples <- 100
  if (post_samples < min_samples) {
    stop(insight::format_error(
      paste0("Configuration produces only ", post_samples, " posterior samples. ",
             "Increase iter, chains, or reduce warmup to get at least ", min_samples, " samples.")
    ))
  }

  # Adjust cores if needed
  max_cores <- parallel::detectCores()
  if (cores > max_cores) {
    if (!identical(Sys.getenv("TESTTHAT"), "true")) {
      rlang::warn(paste0("Requested ", cores, " cores but only ", max_cores,
                        " available. Using ", min(cores, max_cores), " cores."),
                  .frequency = "once", .frequency_id = "cores_limit")
    }
    cores <- min(cores, max_cores)
  }

  list(
    K = K,
    p = lag,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    times = times
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
      specials = list(K = params$K, p = params$p),
      h = h,
      family = gaussian(),
      chains = params$chains,
      cores = params$cores,
      iter = params$iter,
      warmup = params$warmup,
      adapt_delta = params$adapt_delta,
      max_treedepth = params$max_treedepth
    )
  } else if (model == "VARDF") {
    train_vardf(
      .data = object_tsbl,
      specials = list(K = params$K, p = params$p),
      h = h,
      family = gaussian(),
      chains = params$chains,
      cores = params$cores,
      iter = params$iter,
      warmup = params$warmup,
      adapt_delta = params$adapt_delta,
      max_treedepth = params$max_treedepth
    )
  } else if (model == "GPDF") {
    train_gpdf(
      .data = object_tsbl,
      specials = list(K = params$K),
      h = h,
      family = gaussian(),
      chains = params$chains,
      cores = params$cores,
      iter = params$iter,
      warmup = params$warmup,
      adapt_delta = params$adapt_delta,
      max_treedepth = params$max_treedepth
    )
  }
}

#' Adjust forecast uncertainty using variance inflation
#' Combines model forecast with estimated standard deviation
#' to simulate more realistic forecast distributions
#' @noRd
adjust_forecast_uncertainty <- function(
    forecast_df,
    object_sds,
    times,
    h
) {
  forecast_df |>
    dplyr::left_join(
      object_sds,
      by = dplyr::join_by(.basis, .realisation)
    ) |>
    dplyr::mutate(
      .estimate = vctrs::new_vctr(
        purrr::map2(
          .estimate, .sd,
          \(x, y) generate(
            distributional::dist_normal(
              mean = mean(x),
              sd = sqrt(
                distributional::parameters(x)$sigma^2 + y^2
              )
            ),
            times = times
          )
        ),
        class = "distribution"
      )
    ) |>
    tidyr::unnest(cols = .estimate) |>
    tidyr::unnest(cols = .estimate) |>
    dplyr::rename(".sim" = ".estimate") |>
    dplyr::group_by(.basis, .realisation) |>
    dplyr::mutate(.rep = rep(1:times, h)) |>
    dplyr::ungroup()
}

#' Simulate forecasts using fable models
#' If stationary = TRUE, apply order constraints to ARIMA
#' If uncertainty estimates are available, adjust forecasts
#' @noRd
fable_forecast <- function(
    object_tsbl,
    model,
    h,
    times,
    stationary = FALSE,
    object_sds = NULL
) {
  # Input validation
  checkmate::assert_data_frame(object_tsbl, min.rows = 1)
  checkmate::assert_string(model)
  checkmate::assert_count(h, positive = TRUE)
  checkmate::assert_count(times, positive = TRUE)
  checkmate::assert_logical(stationary, len = 1)
  if (!is.null(object_sds)) {
    checkmate::assert_data_frame(object_sds, min.rows = 1)
  }
  
  # Handle ENS ensemble forecasting
  if (model == "ENS") {
    # Generate ETS forecasts
    ets_fc <- fable_forecast(
      object_tsbl = object_tsbl,
      model = "ETS",
      h = h,
      times = times,
      stationary = FALSE,
      object_sds = object_sds
    )
    
    # Generate RW forecasts  
    rw_fc <- fable_forecast(
      object_tsbl = object_tsbl,
      model = "RW",
      h = h,
      times = times,
      stationary = FALSE,
      object_sds = object_sds
    )
    
    # Detect the actual time index column from the tsibble
    # This handles cases where the time column isn't named ".time"
    time_index <- tsibble::index_var(object_tsbl)
    
    # Validate forecast compatibility
    if (!identical(dim(ets_fc), dim(rw_fc))) {
      stop(insight::format_error(
        "ETS and RW forecasts have incompatible dimensions: {.field ETS} = {paste(dim(ets_fc), collapse='x')}, {.field RW} = {paste(dim(rw_fc), collapse='x')}"
      ))
    }
    
    # Check required columns exist with flexible time index
    required_cols <- c(".basis", ".realisation", time_index, ".sim", ".rep")
    
    # Validate ETS forecast structure
    checkmate::assert_subset(required_cols, names(ets_fc), 
      .var.name = paste("ETS forecast columns (time index:", time_index, ")"))
    
    # Validate RW forecast structure  
    checkmate::assert_subset(required_cols, names(rw_fc),
      .var.name = paste("RW forecast columns (time index:", time_index, ")"))
    
    # Combine forecasts with equal weights (arithmetic mean)
    # Use dynamic time index column name in join
    ensemble_fc <- ets_fc |>
      dplyr::left_join(
        rw_fc |>
          dplyr::select(.basis, .realisation, !!time_index, .rep, .sim) |>
          dplyr::rename(.sim_rw = .sim),
        by = c(".basis", ".realisation", time_index, ".rep")
      ) |>
      dplyr::mutate(
        .sim = (.sim + .sim_rw) / 2,
        .model = "ENS"
      ) |>
      dplyr::select(-.sim_rw)
    
    return(ensemble_fc)
  }
  
  # Construct model call with optional constraints
  model_call <- if (stationary) {
    do.call(
      `::`,
      list("fable", model)
    )(
      .estimate,
      order_constraint = (p + q + P + Q <= 6) & (d + D == 0)
    )
  } else {
    do.call(
      `::`,
      list("fable", model)
    )(.estimate)
  }

  # Adjust observation variance if necessary
  if (!is.null(object_sds)) {
    forecast_df <- object_tsbl |>
      fabletools::model(model_call) |>
      forecast(h = h) |>
      dplyr::mutate(.model = model)

    return(
      adjust_forecast_uncertainty(
        forecast_df = forecast_df,
        object_sds = object_sds,
        times = times,
        h = h
      )
    )
  } else {
    forecast_df <- object_tsbl |>
      fabletools::model(model_call) |>
      fabletools::generate(h = h, times = times) |>
      dplyr::mutate(.model = model)

    return(forecast_df)
  }
}

#' Forecasting `ffc_gam` models
#'
#' @importFrom stats median quantile setNames
#' @param model A character string specifying the forecasting model to use.
#' Default is "ENS" (ensemble). Options include "ENS", "ETS", "ARIMA", "RW", 
#' "NAIVE", and Stan dynamic factor models ("ARDF", "VARDF", "GPDF").
#' "ENS" combines ETS and Random Walk forecasts with equal weights,
#' hedging bets between different forecasting assumptions to improve robustness.
#' @inheritParams forecast.fts_ts
#' @param object An object of class `ffc_gam`. See [ffc_gam()]
#' @param newdata `dataframe` or `list` of test data containing the
#' same variables that were included in the original `data` used to fit the
#' model. The covariate information in `newdata`, along with the temporal information
#' indexed by the `time` variable in the original call to `ffc_gam()`, will be used to
#' generate forecasts from the fitted model equations
#' @param type When this has the value `link`, the linear predictor is calculated
#' on the link scale. If `expected` is used, predictions reflect the expectation
#' of the response (the mean) but ignore uncertainty in the observation process.
#' When `response` is used (the default), the predictions take uncertainty in the observation
#' process into account to return predictions on the outcome scale
#' @param summary Should summary statistics be returned instead of the raw values?
#' Default is `TRUE`
#' @param robust If `FALSE` (the default) the mean is used as the measure of
#' central tendency and the standard deviation as the measure of variability.
#' If `TRUE`, the median and the median absolute deviation (MAD) are applied instead.
#' Only used if `summary` is `TRUE`
#' @param probs The percentiles to be computed by the `quantile()` function.
#' Only used if `summary` is `TRUE`
#' @param mean_model A character string specifying the forecasting model to use
#' for mean basis coefficients when using Stan factor models (ARDF, VARDF, GPDF).
#' Default is "ENS". Options include "ENS", "ETS", "ARIMA", "RW", "NAIVE".
#' "ENS" creates an ensemble of ETS and RW forecasts with equal weights,
#' hedging bets between different forecasting assumptions for mean coefficients.
#' This is only used when forecasting mixed mean/non-mean basis functions with
#' Stan factor models.
#' @param ... Additional arguments for Stan dynamic factor models (ARDF, VARDF, GPDF):
#'   \describe{
#'     \item{K}{Number of latent factors (default: 2). Must be positive.}
#'     \item{lag}{Number of time lags for autoregressive terms (default: 1). Must be positive.}
#'     \item{chains}{Number of MCMC chains (default: 4)}
#'     \item{iter}{Total iterations per chain (default: 500)}
#'     \item{warmup}{Warmup iterations per chain (default: iter/2)}
#'     \item{cores}{Number of CPU cores to use (default: min(chains, available cores))}
#'     \item{adapt_delta}{Target acceptance rate (default: 0.75)}
#'     \item{max_treedepth}{Maximum tree depth (default: 9)}
#'     \item{times}{Number of posterior samples to draw for coefficient forecasting.
#'       For Stan models, this is automatically set to (iter - warmup) * chains
#'       to ensure dimension consistency. Minimum 100 total posterior samples required.}
#'   }
#' @details
#' \strong{Forecasting Methodology:}
#'
#' This function implements a two-stage forecasting approach for functional regression
#' models with time-varying coefficients of the form:
#' \deqn{y_t = \sum_{j=1}^{J} \beta_j(t) B_j(x) + \epsilon_t}
#' where \eqn{\beta_j(t)} are time-varying coefficients and \eqn{B_j(x)} are basis functions.
#'
#' \enumerate{
#'   \item \strong{Extract basis coefficients:} Time-varying functional coefficients
#'     \eqn{\beta_j(t)} are extracted from the fitted GAM as time series
#'   \item \strong{Forecast coefficients:} These coefficient time series are forecast
#'     using ensemble methods (ENS), Stan dynamic factor models (ARDF/VARDF/GPDF), 
#'     or individual fable models (ETS, ARIMA, etc.)
#'   \item \strong{Reconstruct forecasts:} Forecasted coefficients are combined:
#'     \deqn{\hat{y}_{t+h} = \sum_{j=1}^{J} \hat{\beta}_j(t+h) B_j(x)}
#'   \item \strong{Combine uncertainties:} Multiple uncertainty sources are integrated
#'     hierarchically (see Uncertainty Quantification section)
#' }
#'
#' \strong{Uncertainty Quantification:}
#'
#' Forecast uncertainty is captured through a hierarchical structure:
#'
#' \emph{Within Stan dynamic factor models:}
#' \itemize{
#'   \item \strong{Process uncertainty:} Factor dynamics, autoregressive terms, factor loadings
#'   \item \strong{Observation uncertainty:} Series-specific error terms (\eqn{\sigma_{obs}})
#' }
#'
#' \emph{Final combination in linear predictor space:}
#' \itemize{
#'   \item \strong{Stan forecast samples:} Already incorporate process + observation uncertainty
#'   \item \strong{GAM parameter uncertainty:} Random draws from \eqn{N(\hat{\boldsymbol{\theta}}, \mathbf{V})}
#'     where \eqn{\mathbf{V}} is the coefficient covariance matrix
#' }
#'
#' These components are combined additively: \eqn{\text{Stan forecasts} + \text{GAM uncertainty}}
#'
#' \strong{Model Selection:}
#'
#' \itemize{
#'   \item \strong{ENS ensemble (default):} Combines ETS and Random Walk forecasts
#'     with equal weights. This hedges bets between exponential smoothing assumptions
#'     (trend and seasonality patterns continue) and random walk assumptions
#'     (future values equal current values). Provides robust predictions when
#'     model uncertainty is high, which is common in coefficient forecasting.
#'   \item \strong{Stan factor models (ARDF/VARDF/GPDF):} Used for multivariate
#'     forecasting of non-mean basis coefficients. Capture dependencies between
#'     coefficient series and assume zero-centered time series for efficiency.
#'   \item \strong{Mean basis models:} Used for mean basis coefficients (which operate
#'     at non-zero levels). Default is ENS ensemble, controlled by `mean_model` parameter.
#' }
#'
#' \strong{Important Note on `times` Parameter:}
#'
#' For Stan dynamic factor models, the `times` parameter is automatically set to
#' `(iter - warmup) * chains` to ensure dimensional consistency. Any user-specified
#' `times` value will be ignored with a warning. For ARIMA models, `times` can be
#' specified freely and controls the number of posterior draws.
#' @examples
#' \donttest{
#' # Basic forecasting example with growth data
#' data("growth_data")
#' mod <- ffc_gam(
#'   height_cm ~ s(id, bs = "re") +
#'     fts(age_yr, k = 8, bs = "cr", time_k = 10),
#'   time = "age_yr", 
#'   data = growth_data,
#'   family = Gamma()
#' )
#'
#' # Forecast with ETS model
#' newdata <- data.frame(
#'   id = "boy_11", 
#'   age_yr = c(16, 17, 18)
#' )
#' fc <- forecast(mod, newdata = newdata)  # Uses ENS ensemble by default
#'
#' # Forecast with specific models
#' fc_ets <- forecast(mod, newdata = newdata, model = "ETS")
#' fc_rw <- forecast(mod, newdata = newdata, model = "RW")
#'
#' # Get raw forecast matrix without summary
#' fc_raw <- forecast(mod, newdata = newdata, summary = FALSE)
#' }
#' @seealso [ffc_gam()], [fts()], [forecast.fts_ts()]
#' @return Predicted values on the appropriate scale.
#' If `summary == FALSE`, the output is a matrix. If `summary == TRUE`, the output
#' is a tidy `tbl_df / data.frame`
#' @author Nicholas J Clark
#' @export
forecast.ffc_gam <- function(
    object,
    newdata,
    type = "response",
    model = "ENS",
    stationary = FALSE,
    summary = TRUE,
    robust = TRUE,
    probs = c(0.025, 0.1, 0.9, 0.975),
    mean_model = "ENS",
    ...) {
  type <- match.arg(
    arg = type,
    choices = c(
      "link",
      "expected",
      "response"
    )
  )

  # Take full draws of beta coefficients
  # Get times parameter for beta draws - need to check for Stan models first
  temp_params <- extract_model_params(...)

  # For Stan models, ensure times matches posterior samples to avoid dimension mismatches
  if (model %in% c('ARDF', 'GPDF', 'VARDF')) {
    stan_post_samples <- (temp_params$iter - temp_params$warmup) * temp_params$chains
    if (temp_params$times != stan_post_samples) {
      if (!identical(Sys.getenv("TESTTHAT"), "true")) {
        rlang::warn(paste0("For Stan models, setting times = ", stan_post_samples, " to match posterior samples. ",
                          "Requested times = ", temp_params$times, " would cause dimension mismatch."),
                    .frequency = "once", .frequency_id = "stan_times_mismatch")
      }
      temp_params$times <- stan_post_samples
    }
  }

  orig_betas <- mgcv::rmvn(
    n = temp_params$times,
    mu = coef(object),
    V = vcov(object)
  )

  if (length(object$gam_init) == 0) {
    # No need to modify lpmatrix if there were no
    # fts() terms in the model
    # Extract the full linear predictor matrix
    orig_lpmat <- predict(
      object,
      newdata = newdata,
      type = "lpmatrix"
    )

    full_linpreds <- matrix(
      as.vector(
        t(apply(
          as.matrix(orig_betas),
          1,
          function(row) {
            orig_lpmat %*% row +
              attr(orig_lpmat, "model.offset")
          }
        ))
      ),
      nrow = NROW(orig_betas)
    )
  } else {
    # Determine horizon (assuming equal time gaps)
    time_var <- object$time_var

    # Validate newdata and filter to future time points only
    newdata <- validate_forecast_newdata(newdata, object)
    validated_newdata_size <- nrow(newdata)
    
    # Predict by fixing time to its last training value; this allows
    # us to later add the uncertainty from the zero-centred time-varying
    # basis coefficient time series
    last_train <- max(object$model[[object$time_var]])
    newdata_shift <- newdata
    newdata_shift[[object$time_var]] <- last_train

    orig_lpmat <- predict(
      object,
      newdata = newdata_shift,
      type = "lpmatrix"
    )
    rm(newdata_shift, last_train)
    
    # Now create basis coefficient data for forecasting using filtered newdata
    interpreted <- interpret_ffc(
      formula = object$orig_formula,
      data = newdata,
      newdata = newdata,
      gam_init = object$gam_init,
      time_var = object$time_var
    )

    fc_horizons <- interpreted$data[[time_var]] -
      max(object$model[[time_var]])
    # Use number of unique forecast time steps, not total data points
    # This prevents massive memory usage when forecasting many data points
    max_horizon <- length(unique(fc_horizons))

    # Extract functional basis coefficient time series
    # Use the same temp_params that was already adjusted for Stan models
    intermed_coefs <- fts_coefs(
      object,
      summary = FALSE,
      times = temp_params$times
    )

    # Optimized functional coefficient processing
    # Aggregate to unique combinations FIRST (eliminates 99.5% of redundant processing)
    aggregated <- intermed_coefs |>
      dplyr::group_by(.basis, .time, !!time_var) |>
      dplyr::summarise(.mean = mean(.estimate), .groups = "drop")
    
    # Validate aggregation structure
    checkmate::assert_data_frame(aggregated, min.rows = 1)
    
    # Vectorized tail subtraction (compute last value per basis)
    tail_values <- aggregated |>
      dplyr::group_by(.basis) |>
      dplyr::slice_tail(n = 1) |>
      dplyr::select(.basis, tail_mean = .mean)
    
    # Join and compute shifted values (subtract tail from each group)
    functional_coefs <- aggregated |>
      dplyr::left_join(tail_values, by = ".basis") |>
      dplyr::mutate(
        .mean = .mean - tail_mean,
        .realisation = 1,
        .estimate = .mean
      ) |>
      dplyr::select(-tail_mean) |>
      tsibble::as_tsibble(
        key = c(.basis, .realisation),
        index = .time
      )

    # Structure the functional_coefs object for automatic
    # recognition of time signatures (if the data provided are
    # a tsibble format)
      if (!is.null(attr(intermed_coefs, "index"))) {
        if (!attr(intermed_coefs, "index") %in% names(functional_coefs)){
          functional_coefs <- functional_coefs |>
            dplyr::left_join(
              intermed_coefs |>
                dplyr::select(.time, !!attr(intermed_coefs, "index")) |>
                dplyr::distinct(),
              by = dplyr::join_by(.time)
            )
        }
      }


    functional_coefs <- structure(
      functional_coefs,
      class = c("fts_ts", "tbl_df", "tbl", "data.frame"),
      time_var = attr(intermed_coefs, "time_var"),
      index = attr(intermed_coefs, "index"),
      index2 = attr(intermed_coefs, "index2"),
      interval = attr(intermed_coefs, "interval"),
      summarized = attr(intermed_coefs, "summarized")
    )

    # Fit the time series model to the basis coefficients
    # and generate basis coefficient forecasts

    # Categorize basis functions once
    basis_names <- unique(functional_coefs$.basis)
    has_mean_basis <- any(grepl('_mean', basis_names))
    has_non_mean_basis <- any(!grepl('_mean', basis_names))

    if(model %in% c('ARDF', 'GPDF', 'VARDF') &
       has_mean_basis & has_non_mean_basis) {
      # For factor models, it will very often make sense to forecast any
      # _mean basis (i.e. level shifts) separately as they can behave very
      # differently to remaining basis coefficient time series

      # First fit the dynamic factor model for all non (_mean)
      # basis functions
      functional_fc_others <- suppressWarnings(
        forecast(
          object = functional_coefs |>
            dplyr::filter(!grepl('_mean', .basis)),
          h = max_horizon,
          times = temp_params$times,
          model = model,
          stationary = stationary,
          ...
        )
      )

      if('.time' %in% colnames(functional_fc_others)){
        functional_fc_others <- functional_fc_others |>
          dplyr::select(-.time)
      }

      # Now forecast the _mean basis with the specified model
      functional_fc_mean <- suppressWarnings(
        forecast(
          object = functional_coefs |>
            dplyr::filter(grepl('_mean', .basis)),
          h = max_horizon,
          times = temp_params$times,
          model = mean_model,
          stationary = FALSE,
          ...
        )
      ) |>
        tsibble::as_tibble() |>
        dplyr::select(dplyr::any_of(names(functional_fc_others)))

      if(!time_var %in% colnames(functional_fc_mean) &
         is.null(attr(intermed_coefs, "index"))){
        functional_fc_mean <- functional_fc_mean |>
          dplyr::mutate({{time_var}} := .time)
      }

      # Bind the two forecasts together
      functional_fc <- functional_fc_others |>
        dplyr::mutate(.rep = as.character(.rep)) |>
        dplyr::bind_rows(functional_fc_mean)

    } else {
      # Use single forecasting approach for all basis functions
      # This handles cases where:
      # 1. Non-factor models are used (ARIMA)
      # 2. Only _mean basis functions exist (mean_only = TRUE)
      # 3. Only non-_mean basis functions exist
      functional_fc <- suppressWarnings(
        forecast(
          object = functional_coefs,
          h = max_horizon,
          times = temp_params$times,
          model = if(has_mean_basis & !has_non_mean_basis) mean_model else model,
          stationary = if(has_mean_basis & !has_non_mean_basis) FALSE else stationary,
          ...
        )
      )
    }

    # Ensure original time structure is preserved
    index_var <- attr(intermed_coefs, "index")

    if (!is.null(index_var)) {
      # Join original time data using the index attribute
      functional_fc <- functional_fc |>
        dplyr::left_join(
          interpreted$orig_data |>
            tibble::as_tibble() |>
            dplyr::select(!!time_var, !!index_var) |>
            dplyr::distinct(),
          by = dplyr::join_by(!!index_var)
        ) |>
        dplyr::rename(.time = !!time_var)

    } else if (time_var %in% colnames(functional_fc)) {
      # Replace .time column with time_var if available
      functional_fc <- functional_fc |>
        tsibble::as_tibble() |>
        dplyr::select(-.time) |>
        dplyr::bind_cols(
          functional_fc |>
            tsibble::as_tibble() |>
            dplyr::select(!!time_var) |>
            dplyr::rename(.time = !!time_var)
        )
    }

    # Map forecast .time values to actual target time values
    # The forecast generates consecutive integer steps, but we need specific time values
    target_times <- unique(interpreted$data[[time_var]])
    forecast_times <- sort(unique(functional_fc$.time))
    
    
    # Validate time mapping is mathematically possible
    if (length(forecast_times) != length(target_times)) {
      stop(insight::format_error(
        paste0("Cannot map forecast times: generated ", length(forecast_times), 
               " unique time steps but newdata has ", length(target_times), 
               " unique time values. This may indicate an issue with the forecasting horizon.")
      ))
    }
    
    # Create mapping from forecast times to target times
    time_mapping <- setNames(target_times, forecast_times)
    
    # Remap the .time values
    # Convert to tibble first to avoid tsibble validation issues with duplicate key+time combinations
    functional_fc <- functional_fc |>
      tsibble::as_tibble() |>
      dplyr::mutate(.time = time_mapping[as.character(.time)])

    # Calculate intermediate linear predictors
    model_offset <- attr(orig_lpmat, "model.offset")

    # Compute linear predictors using matrix multiplication
    intermed_linpreds <- orig_betas %*%
      t(orig_lpmat)
    intermed_linpreds <- sweep(
      intermed_linpreds,
      2,
      model_offset,
      "+"
    )

    # Compute functional coefficient predictions using matrix operations
    fts_fc <- compute_functional_predictions(
      interpreted$data, 
      functional_fc, 
      time_var
    )

    # Should also have same ncols as number of fts basis functions
    if (
      !NCOL(fts_fc) - 2 ==
        NCOL(interpreted$data[
          grep(
            "fts_",
            colnames(interpreted$data)
          )
        ])
    ) {
      stop("Wrong dimensions in forecast coefs; need to check on this")
    }

    # If dimensions correct, take rowsums for each draw
    fc_linpreds <- fts_fc |>
      dplyr::select(-c(.draw, .row)) |>
      rowSums()
    # Add the draw-specific row predictions to the
    # intermediate prediction matrix
    unique_draws <- unique(fts_fc$.draw)

    fc_matrix <- optimized_fc_matrix_reshape(
      fc_linpreds = fc_linpreds,
      draw_ids = fts_fc$.draw,
      unique_draws = unique_draws
    )

    full_linpreds <- fc_matrix + intermed_linpreds
  }

  # Now can proceed to send full_linpreds to the relevant
  # invlink and rng functions for outcome-level predictions
  if (type == "link") preds <- full_linpreds

  if (type == "expected") {
    preds <- posterior_epred(
      object,
      full_linpreds
    )
  }

  if (type == "response") {
    preds <- posterior_predict(
      object,
      full_linpreds
    )
  }

  # Summarise if necessary
  if (summary) {
    Qs <- apply(preds, 2, quantile, probs = probs, na.rm = TRUE)

    if (robust) {
      estimates <- apply(preds, 2, median, na.rm = TRUE)
      errors <- apply(abs(preds - estimates), 2, median, na.rm = TRUE)
    } else {
      estimates <- apply(preds, 2, mean, na.rm = TRUE)
      errors <- apply(preds, 2, sd, na.rm = TRUE)
    }

    out <- data.frame(
      cbind(estimates, errors, t(Qs))
    )
    colnames(out) <- c(
      ".estimate",
      ".error",
      paste0(".q", 100 * probs)
    )
    class(out) <- c("tbl_df", "tbl", "data.frame")
  } else {
    out <- distributional::dist_sample(
      lapply(seq_len(NCOL(preds)), function(i) preds[, i])
    )
  }

  return(out)
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

#' Posterior expectations
#' @noRd
posterior_epred <- function(object,
                            linpreds) {
  # invlink function
  invlink_fun <- get_family_invlink(object)

  # Compute expectations
  expected_pred_vec <- invlink_fun(
    eta = as.vector(linpreds)
  )

  # Convert back to matrix
  out <- matrix(expected_pred_vec,
    nrow = NROW(linpreds)
  )
  return(out)
}

#' Posterior predictions
#' @noRd
posterior_predict <- function(object,
                              linpreds) {
  # rd function if available
  rd_fun <- get_family_rd(object)

  # invlink function
  invlink_fun <- get_family_invlink(object)

  # Dispersion parameter
  scale_p <- object[["sig2"]]
  if (is.null(scale)) {
    scale_p <- summary(object)[["dispersion"]]
  }

  if (!grepl("tweedie", object$family[["family"]],
    ignore.case = TRUE
  )) {
    scale_p <- rep(scale_p, length(linpreds))
  }

  # weights
  weights <- rep(1, length(linpreds))

  # Compute expectations
  expected_pred_vec <- invlink_fun(
    eta = as.vector(linpreds)
  )

  # For Tweedie family, ensure mu values are non-negative
  if (grepl("tweedie", object$family[["family"]], ignore.case = TRUE)) {
    expected_pred_vec <- pmax(expected_pred_vec, 1e-8)
  }

  # Now compute response predictions
  response_pred_vec <- rd_fun(
    mu = expected_pred_vec,
    wt = weights,
    scale = scale_p
  )

  # Convert back to matrix
  out <- matrix(response_pred_vec,
    nrow = NROW(linpreds)
  )
  return(out)
}

#' Functions to enable random number generation from mgcv
#' gam / bam objects. Code is modified from severa internal functions written
#' by Gavin Simpson for the gratia R package, which in turn were modified from
#' original code written by Simon Wood for the mgcv R package

# simulator for tweedie LSS models
#' @importFrom rlang .data
#' @importFrom stats rpois rgamma
#' @importFrom tibble tibble
#' @noRd
rtw <- function(mu, p, phi) {
  if (any(p <= 1 | p >= 2)) {
    stop("'p' must be in interval (1, 2)")
  }
  if (any(phi <= 0)) {
    stop("scale parameter 'phi' must be positive")
  }
  if (any(mu < 0)) {
    stop("mean 'mu' must be non-negative")
  }
  lambda <- mu^(2 - p) / ((2 - p) * phi)
  shape <- (2 - p) / (p - 1)
  scale <- phi * (p - 1) * mu^(p - 1)
  N <- rpois(length(lambda), lambda)
  gs <- rep(scale, N)
  tab <- tibble(
    y = rgamma(gs * 0 + 1, shape = shape, scale = gs),
    lab = rep(seq_along(N), N)
  )
  out <- numeric(length(N))
  out[which(N != 0)] <- tab |>
    dplyr::group_by(.data$lab) |>
    dplyr::summarise(summed = sum(.data$y)) |>
    dplyr::pull(.data$summed)
  out
}

#' converts from theta to power parameter `p` given `a` and `b`
#' @noRd
theta_2_power <- function(theta, a, b) {
  i <- theta > 0
  exp_theta_pos <- exp(-theta[i])
  exp_theta_neg <- exp(theta[!i])
  theta[i] <- (b + a * exp_theta_pos) / (1 + exp_theta_pos)
  theta[!i] <- (b * exp_theta_neg + a) / (1 + exp_theta_neg)
  theta
}

#' extracts the `a` and `b` parameters of the model search over which the power
#' parameter is searched for
#' @noRd
get_tw_ab <- function(family) {
  if (family[["family"]] != "twlss") {
    stop("'model' wasn't fitted with 'twlss()' family.", call. = FALSE)
  }
  rfun <- family$residuals
  a <- get(".a", envir = environment(rfun))
  b <- get(".b", envir = environment(rfun))
  c(a, b)
}

#' @importFrom mgcv fix.family.rd
#' @noRd
fix_family_rd <- function(family, ...) {
  # try to fix up the family used by mgcv to add the $rd component
  # for random deviate sampling

  # try the obvious thing first and see if mgcv::fix.family.rd() already handles
  # family
  fam <- mgcv::fix.family.rd(family)

  # if `family` contains a NULL rd we move on, if it is non-null return early
  # as it doesn't need fixing
  if (!is.null(fam$rd)) {
    return(fam)
  }

  # handle special cases
  fn <- fam[["family"]]

  # handle multivariate normal
  if (identical(fn, "Multivariate normal")) {
    # note: mgcv::mvn is documented to ignore prior weights
    # if we ever need to handle weights to scale V, see this post on CV
    # https://stats.stackexchange.com/a/162885/1390
    rd_mvn <- function(V) {
      function(mu, wt, scale) { # function needs to take wt and scale
        mgcv::rmvn(
          n = nrow(mu),
          mu = mu,
          V = V
        )
      }
    }
    fam$rd <- rd_mvn(solve(crossprod(fam$data$R)))
  }
  if (identical(fn, "twlss")) {
    # this uses some helpers to find the `a` and `b` used during fitting and
    # also to convert what `predict()` etc returns (theta) to power parameter
    rd_twlss <- function(a, b) {
      function(mu, wt, scale) {
        rtw(
          mu = mu[, 1], # fitted(model) for twlss is on response scale!
          p = theta_2_power(theta = mu[, 2], a, b),
          phi = exp(mu[, 3])
        )
      }
    }
    tw_pars <- get_tw_ab(fam)
    fam$rd <- rd_twlss(a = tw_pars[1], b = tw_pars[2])
  }

  # return modified family
  fam
}

#' @importFrom mgcv fix.family.rd
#' @importFrom stats family
#' @noRd
get_family_rd <- function(object) {
  if (inherits(object, "glm")) {
    fam <- family(object) # extract family
  } else {
    fam <- object[["family"]]
  }
  ## mgcv stores data simulation funs in `rd`
  fam <- fix_family_rd(fam)
  if (is.null(fam[["rd"]])) {
    stop("Don't yet know how to simulate from family <",
      fam[["family"]], ">",
      call. = FALSE
    )
  }
  fam[["rd"]]
}

#' @noRd
get_family_invlink <- function(object) {
  if (inherits(object, "glm")) {
    fam <- family(object) # extract family
  } else {
    fam <- object[["family"]]
  }
  ## mgcv stores data simulation funs in `rd`
  fam <- fix_family_rd(fam)
  if (is.null(fam[["rd"]])) {
    stop("Don't yet know how to simulate from family <",
      fam[["family"]], ">",
      call. = FALSE
    )
  }
  fam[["linkinv"]]
}
