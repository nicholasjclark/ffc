#' Forecasting functional basis coefficients (Internal)
#'
#' @description
#' \strong{For internal use:} This function is primarily used internally by
#' `forecast.ffc_gam()`. Most users should call `forecast()` directly on their
#' `ffc_gam` object instead of using this function, which requires properly
#' structured coefficient data and has specific format requirements.
#'
#' @param object An object of class `fts_ts` containing time-varying
#' basis function coefficients extracted from an `ffc_gam` object using [fts_coefs()]
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
  silent <- if ("silent" %in% names(dots)) dots$silent else get_stan_param("silent")
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
  checkmate::assert_logical(silent, len = 1)
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
    silent = silent,
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
      max_treedepth = params$max_treedepth,
      silent = params$silent
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
      max_treedepth = params$max_treedepth,
      silent = params$silent
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
      max_treedepth = params$max_treedepth,
      silent = params$silent
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
        "ETS and RW forecasts have incompatible dimensions: ETS = {paste(dim(ets_fc), collapse='x')}, RW = {paste(dim(rw_fc), collapse='x')}"
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
#' @param ... Additional arguments for Stan dynamic factor models (ARDF, VARDF, GPDF).
#'   Key arguments include: K (number of factors, default: 2), lag (AR order, default: 1),
#'   chains (MCMC chains, default: 4), iter (iterations, default: 500), 
#'   silent (suppress progress, default: TRUE), cores, adapt_delta, max_treedepth.
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
  
  # Store original row order before any processing
  # This will be used to restore the original ordering at the end
  original_newdata_rows <- nrow(newdata)
  if (!".original_row_id" %in% names(newdata)) {
    newdata$.original_row_id <- seq_len(original_newdata_rows)
  }
  
  # Initialize surviving_rows for row order restoration
  surviving_rows <- NULL

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
    # No fts() terms - still need to track row IDs for restoration
    surviving_rows <- newdata$.original_row_id
    
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
    # This may filter out rows and reorder data
    newdata_validated <- validate_forecast_newdata(newdata, object)
    validated_newdata_size <- nrow(newdata_validated)
    
    # Extract which original rows survived validation
    if (".original_row_id" %in% names(newdata_validated)) {
      surviving_rows <- newdata_validated$.original_row_id
    } else {
      # Fallback if validation didn't preserve IDs
      surviving_rows <- NULL
    }
    
    # Use validated data for forecasting
    newdata <- newdata_validated
    
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
        )
      
      # Handle duplicate .time column
      if (".time" %in% colnames(functional_fc)) {
        # Remove duplicate time_var column since .time already exists
        functional_fc <- functional_fc |>
          dplyr::select(-!!time_var)
      } else {
        # Rename time_var to .time if .time doesn't exist
        functional_fc <- functional_fc |>
          dplyr::rename(.time = !!time_var)
      }

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

    # Map forecast times to target times based on sorted temporal order
    # This preserves relative temporal ordering while mapping to target values
    target_times <- interpreted$data |>
      dplyr::distinct(!!rlang::sym(time_var)) |>
      dplyr::arrange(!!rlang::sym(time_var)) |>
      dplyr::pull(!!rlang::sym(time_var))
    
    forecast_times <- functional_fc |>
      dplyr::distinct(.time) |>
      dplyr::arrange(.time) |>
      dplyr::pull(.time)
    
    # Floating point tolerance for numeric time comparisons
    TIME_TOLERANCE <- 1e-10
    
    # Create time mapping based on sorted temporal order
    if (length(forecast_times) == length(target_times)) {
      # Direct 1:1 mapping when counts match
      time_mapping <- data.frame(
        .time = forecast_times,
        .time_target = target_times
      )
    } else if (length(forecast_times) < length(target_times)) {
      # Forecast times are subset of target times - map by temporal value
      if (is.numeric(forecast_times) && is.numeric(target_times)) {
        # Use tolerance for floating point comparison
        available_targets <- target_times[sapply(target_times, function(t) 
          any(abs(forecast_times - t) < TIME_TOLERANCE))]
      } else {
        available_targets <- target_times[target_times %in% forecast_times]
      }
      
      if (length(available_targets) == length(forecast_times)) {
        time_mapping <- data.frame(
          .time = forecast_times,
          .time_target = forecast_times  # Direct mapping for matching values
        )
      } else {
        stop(insight::format_error(
          paste0("Cannot map forecast times to target times. ",
                 "Forecast generated {.field ", length(forecast_times), 
                 "} time steps but only {.field ", length(available_targets),
                 "} match target times.")
        ), call. = FALSE)
      }
    } else {
      # More forecast times than target times - internal error
      stop(insight::format_error(
        paste0("Generated {.field ", length(forecast_times), 
               "} forecast times but found {.field ", length(target_times), 
               "} target times. This suggests an internal forecasting error.")
      ), call. = FALSE)
    }
    
    # Remap the .time values using the mapping table
    functional_fc <- functional_fc |>
      tsibble::as_tibble() |>
      dplyr::left_join(time_mapping, by = ".time") |>
      dplyr::mutate(.time = .time_target) |>
      dplyr::select(-.time_target)

    # Calculate intermediate linear predictors  
    model_offset <- attr(orig_lpmat, "model.offset")

    # Handle multi-parameter offset expansion for distributional models
    # Based on actual mgcv offset structure discovered through debugging
    if (is.list(model_offset)) {
      # Multi-parameter case: expand offsets using lpi from lpmatrix attributes
      lpi <- attr(orig_lpmat, "lpi")
      if (is.null(lpi)) {
        stop(insight::format_error(
          "Multi-parameter model missing {.field lpi} attribute in lpmatrix for offset expansion"
        ), call. = FALSE)
      }
      
      checkmate::assert_list(model_offset, types = "numeric", min.len = 1)
      checkmate::assert_list(lpi, types = "numeric", min.len = 1)
      
      # Validate that offset list length matches lpi list length
      if (length(model_offset) != length(lpi)) {
        stop(insight::format_error(
          paste0("Offset list length {.field ", length(model_offset), 
                 "} doesn't match parameter count {.field ", length(lpi), "}")
        ), call. = FALSE)
      }
      
      expanded_offset <- numeric(ncol(orig_lpmat))
      
      # Expand each parameter's offset to its coefficient indices
      for (i in seq_along(model_offset)) {
        if (!is.null(model_offset[[i]]) && length(model_offset[[i]]) > 0) {
          # Use the coefficient indices from lpi for this parameter
          coef_indices <- lpi[[i]]
          expanded_offset[coef_indices] <- model_offset[[i]]
        }
      }
      model_offset <- expanded_offset
      
    } else if (is.null(model_offset)) {
      model_offset <- rep(0, nrow(orig_lpmat))
    } else {
      # Single parameter case - handle scalar or per-observation offsets
      checkmate::assert_numeric(model_offset)
      if (length(model_offset) == 1) {
        # Scalar offset - replicate for all observations
        model_offset <- rep(model_offset, nrow(orig_lpmat))
      } else if (length(model_offset) == nrow(orig_lpmat)) {
        # Already correct length (per-observation) - use as is
        # No modification needed
      } else {
        # Unexpected length - use first element as scalar fallback
        model_offset <- rep(model_offset[1], nrow(orig_lpmat))
      }
    }

    # Compute linear predictors using matrix multiplication
    # For distributional families, compute parameter-specific linear predictors
    lpi_attr <- attr(orig_lpmat, "lpi")
    
    if (!is.null(lpi_attr) && is_distributional_family(object$family)) {
      # Distributional family: compute parameter-specific linear predictors
      
      # Apply offset if needed
      lpmat_to_use <- orig_lpmat
      if (!is.null(model_offset) && !all(model_offset == 0)) {
        for (i in seq_len(nrow(orig_lpmat))) {
          lpmat_to_use[i, ] <- orig_lpmat[i, ] + model_offset
        }
      }
      
      # Compute parameter-specific linear predictors
      n_params <- object$family$nlp
      n_time_points <- nrow(lpmat_to_use)
      n_draws <- nrow(orig_betas)
      
      # Result matrix: rows = draws, cols = (n_params * n_time_points)
      intermed_linpreds <- matrix(0, nrow = n_draws, ncol = n_params * n_time_points)
      
      for (p in seq_len(n_params)) {
        param_coef_indices <- lpi_attr[[p]]
        
        # Extract parameter-specific coefficients and design matrix
        param_betas <- orig_betas[, param_coef_indices, drop = FALSE]
        param_lpmat <- lpmat_to_use[, param_coef_indices, drop = FALSE]
        
        # Compute linear predictors for this parameter  
        param_linpreds <- param_betas %*% t(param_lpmat)
        
        # Store in appropriate columns
        col_start <- (p - 1) * n_time_points + 1
        col_end <- p * n_time_points
        intermed_linpreds[, col_start:col_end] <- param_linpreds
      }
      
    } else {
      # Single parameter family: use existing computation
      if (is.null(model_offset) || all(model_offset == 0)) {
        # No offset needed
        intermed_linpreds <- orig_betas %*% t(orig_lpmat)
      } else {
        # Apply offset by adding to linear predictor matrix before multiplication
        # Each row of orig_lpmat gets the offset added
        orig_lpmat_with_offset <- orig_lpmat
        for (i in seq_len(nrow(orig_lpmat))) {
          orig_lpmat_with_offset[i, ] <- orig_lpmat[i, ] + model_offset
        }
        intermed_linpreds <- orig_betas %*% t(orig_lpmat_with_offset)
      }
    }

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

    # If dimensions correct, compute parameter-specific functional predictions
    # For distributional families, need to separate by parameter
    original_row_indices <- fts_fc$.row
    unique_draws <- unique(fts_fc$.draw)
    
    # Input validation
    checkmate::assert_data_frame(fts_fc, min.rows = 1)
    checkmate::assert_true(all(c(".draw", ".row") %in% names(fts_fc)))
    
    if (is_distributional_family(object$family)) {
      # Parameter-specific functional coefficient computation
      n_params <- object$family$nlp
      n_time_points <- length(unique(original_row_indices))
      n_draws <- length(unique_draws)
      
      # Validation
      checkmate::assert_count(n_params, positive = TRUE)
      checkmate::assert_count(n_time_points, positive = TRUE)
      checkmate::assert_count(n_draws, positive = TRUE)
      
      # Initialize parameter-specific fc_matrix: draws × (params * time_points)
      fc_matrix <- matrix(0, nrow = n_draws, ncol = n_params * n_time_points)
      
      # Parameter names for column identification
      param_names <- c("location", "scale", "shape")
      
      for (p in seq_len(n_params)) {
        param_name <- param_names[p]
        
        # Select columns for this parameter using established pattern
        param_cols <- grep(paste0("^", param_name, "_fts_"), names(fts_fc), 
                          value = TRUE)
        
        # Validation: ensure parameter columns exist
        if (length(param_cols) == 0) {
          stop(insight::format_error(
            paste0("No functional coefficient columns found for parameter ", 
                   "{.field ", param_name, "}. Expected columns matching ", 
                   "pattern {.field ", param_name, "_fts_*}")
          ), call. = FALSE)
        }
        
        # Compute rowsums for this parameter only
        param_fc_data <- fts_fc |>
          dplyr::select(dplyr::all_of(c(param_cols, ".draw", ".row")))
        
        param_fc_linpreds <- param_fc_data |>
          dplyr::select(-c(.draw, .row)) |>
          rowSums()
        
        # Reshape to matrix form for this parameter
        param_fc_matrix <- optimized_fc_matrix_reshape(
          fc_linpreds = param_fc_linpreds,
          draw_ids = param_fc_data$.draw,
          unique_draws = unique_draws
        )
        
        # Validation: ensure matrix dimensions are correct
        checkmate::assert_matrix(param_fc_matrix, nrows = n_draws, 
                                ncols = n_time_points)
        
        # Store in layout: [param1_times, param2_times, ...]
        col_start <- (p - 1) * n_time_points + 1
        col_end <- p * n_time_points
        fc_matrix[, col_start:col_end] <- param_fc_matrix
      }
    } else {
      # Single parameter: use existing logic unchanged
      fc_linpreds <- fts_fc |>
        dplyr::select(-c(.draw, .row)) |>
        rowSums()
      
      fc_matrix <- optimized_fc_matrix_reshape(
        fc_linpreds = fc_linpreds,
        draw_ids = fts_fc$.draw,
        unique_draws = unique_draws
      )
    }

    
    full_linpreds <- fc_matrix + intermed_linpreds
    
    # Preserve lpi attribute for distributional families
    # For distributional families, the lpi structure needs to match the new matrix layout
    if (is_distributional_family(object$family) && !is.null(attr(orig_lpmat, "lpi"))) {
      # Reconstruct lpi for the new matrix structure: [param1_times, param2_times, ...]
      n_params <- object$family$nlp
      n_time_points <- ncol(full_linpreds) %/% n_params
      
      forecast_lpi <- vector("list", n_params)
      for (p in seq_len(n_params)) {
        col_start <- (p - 1) * n_time_points + 1
        col_end <- p * n_time_points
        forecast_lpi[[p]] <- col_start:col_end
      }
      attr(full_linpreds, "lpi") <- forecast_lpi
    } else if (!is.null(attr(orig_lpmat, "lpi"))) {
      attr(full_linpreds, "lpi") <- attr(orig_lpmat, "lpi")
    }
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
    # Create distribution object with one distribution per observation (row)
    # Each distribution contains all draws for that observation
    out <- distributional::dist_sample(
      lapply(seq_len(NROW(preds)), function(i) preds[i, ])
    )
  }
  
  # Restore original row order if we have tracking information
  # Check that surviving_rows exists and has the right length
  if (!is.null(surviving_rows)) {
    if (length(surviving_rows) == NROW(out)) {
      # Create mapping to restore original order
      row_order <- order(surviving_rows)
      
      if (summary) {
        out <- out[row_order, , drop = FALSE]
      } else {
        # For distributional objects, need to reorder the list
        out <- out[row_order]
      }
    } else {
      # Log mismatch for debugging
      if (!identical(Sys.getenv("TESTTHAT"), "true")) {
        rlang::warn(paste0("Row count mismatch in order restoration: ",
                          "expected ", NROW(out), " rows but found ",
                          length(surviving_rows), " IDs"),
                    .frequency = "once", .frequency_id = "row_count_mismatch")
      }
    }
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

#' Check if family supports multiple parameters
#' @param family mgcv family object
#' @return Logical indicating if family has multiple parameters
#' @noRd
is_distributional_family <- function(family) {
  checkmate::assert_class(family, "family")
  !is.null(family$nlp) && family$nlp > 1
}

#' Extract parameter information from linear predictor matrix
#' @param linpreds Linear predictor matrix from predict() call
#' @param family mgcv family object
#' @return List with parameter information and indices
#' @noRd  
extract_parameter_info_from_lpmat <- function(linpreds, family) {
  checkmate::assert_matrix(linpreds)
  checkmate::assert_class(family, "family")
  
  if (is_distributional_family(family)) {
    # Get parameter indices from lpmatrix attributes (correct mgcv pattern)
    lpi <- attr(linpreds, "lpi")
    if (is.null(lpi)) {
      stop(insight::format_error(
        "Multi-parameter model missing {.field lpi} attribute in lpmatrix"
      ), call. = FALSE)
    }
    
    checkmate::assert_list(lpi, types = "numeric", min.len = 1)
    checkmate::assert_true(length(lpi) == family$nlp)
    
    # Create parameter names based on standard distributional conventions
    par_names <- c("location", "scale", "shape")[seq_len(family$nlp)]
    if (family$nlp > 3) {
      par_names <- c(par_names, paste0("parameter_", 4:family$nlp))
    }
    
    list(
      n_parameters = family$nlp,
      parameter_names = par_names,
      parameter_indices = lpi
    )
  } else {
    list(
      n_parameters = 1L, 
      parameter_names = "location", 
      parameter_indices = NULL
    )
  }
}

#' Split linear predictors by parameter using lpi indices
#' @param linpreds Linear predictor matrix
#' @param parameter_info Parameter information from extract_parameter_info
#' @return List of parameter-specific linear predictors
#' @noRd
split_linear_predictors_by_lpi <- function(linpreds, parameter_info) {
  checkmate::assert_matrix(linpreds)
  checkmate::assert_list(parameter_info, names = "strict")
  checkmate::assert_names(names(parameter_info), 
                         must.include = c("n_parameters", "parameter_indices"))
  
  if (parameter_info$n_parameters == 1L) {
    return(list(location = linpreds))
  }
  
  # Multi-parameter: split using lpi indices
  par_predictions <- vector("list", parameter_info$n_parameters)
  for (i in seq_len(parameter_info$n_parameters)) {
    indices <- parameter_info$parameter_indices[[i]]
    checkmate::assert_integerish(indices, lower = 1, 
                                upper = ncol(linpreds))
    par_predictions[[i]] <- linpreds[, indices, drop = FALSE]
  }
  names(par_predictions) <- parameter_info$parameter_names
  return(par_predictions)
}

#' Apply inverse link functions for distributional families
#' @param par_predictions List of parameter-specific linear predictors
#' @param family mgcv family object
#' @return List of fitted values for each parameter
#' @noRd
apply_distributional_inverse_links <- function(par_predictions, family) {
  checkmate::assert_list(par_predictions, min.len = 1)
  checkmate::assert_class(family, "family")
  
  fitted_pars <- vector("list", length(par_predictions))
  
  if (!is.null(family$linfo)) {
    # Distributional families - use parameter-specific inverse links
    checkmate::assert_true(
      length(par_predictions) <= family$nlp,
      .var.name = "par_predictions length vs family parameters"
    )
    
    for (i in seq_along(par_predictions)) {
      # Defensive check for parameter bounds
      checkmate::assert_true(i <= length(family$linfo))
      fitted_pars[[i]] <- family$linfo[[i]]$linkinv(par_predictions[[i]])
    }
  } else {
    # Standard families - use single inverse link
    checkmate::assert_true(!is.null(family$linkinv))
    
    for (i in seq_along(par_predictions)) {
      fitted_pars[[i]] <- family$linkinv(par_predictions[[i]])
    }
  }
  
  names(fitted_pars) <- names(par_predictions)
  return(fitted_pars)
}

#' Posterior expectations with distributional family support
#' @param object Fitted model object  
#' @param linpreds Linear predictor matrix from predict() call
#' @return Matrix of expectations (location parameter only)
#' @noRd
posterior_epred <- function(object, linpreds) {
  checkmate::assert_matrix(linpreds)
  family <- object$family
  checkmate::assert_class(family, "family")
  
  # For distributional families, expectation = location parameter only
  # gaulss: E[Y] = μ, twlss: E[Y] = μ, betar: E[Y] = μ
  if (is_distributional_family(family)) {
    # Extract parameter information using correct mgcv pattern
    parameter_info <- extract_parameter_info_from_lpmat(linpreds, family)
    
    # Get location parameter (first parameter) indices
    location_indices <- parameter_info$parameter_indices[[1]]
    checkmate::assert_integerish(location_indices, lower = 1, 
                                upper = ncol(linpreds))
    
    # Extract location linear predictors
    location_linpreds <- linpreds[, location_indices, drop = FALSE]
    
    # Apply inverse link to location parameter only
    expectations <- family$linkinv(location_linpreds)
    
  } else {
    # Single parameter family - standard approach
    expectations <- family$linkinv(linpreds)
  }
  
  # Return as matrix
  return(expectations)
}

#' Posterior predictions with optional parameter matrices for distributional families
#' @param object Fitted model object
#' @param linpreds Linear predictor matrix from predict() call
#' @param location_matrix Optional matrix of location parameter values (NULL = extract from linpreds)
#' @param scale_matrix Optional matrix of scale parameter values (NULL = extract from linpreds)
#' @param shape_matrix Optional matrix of shape parameter values (NULL = extract from linpreds)
#' @return Matrix of posterior predictions with nrow = nrow(linpreds)
#' @details For distributional families, parameter matrices allow direct specification
#'   of parameter values, bypassing lpi extraction from linpreds. All parameter matrices
#'   must have same number of rows and columns. When NULL, existing parameter extraction
#'   logic is used. Parameter matrices are validated for correct dimensions.
#' @noRd
posterior_predict <- function(object, linpreds, location_matrix = NULL,
                             scale_matrix = NULL, shape_matrix = NULL) {
  checkmate::assert_matrix(linpreds)
  family <- object$family
  checkmate::assert_class(family, "family")
  
  # Validate parameter matrices if provided
  if (!is.null(location_matrix)) {
    checkmate::assert_matrix(location_matrix, nrows = nrow(linpreds))
  }
  if (!is.null(scale_matrix)) {
    checkmate::assert_matrix(scale_matrix, nrows = nrow(linpreds))
  }
  if (!is.null(shape_matrix)) {
    checkmate::assert_matrix(shape_matrix, nrows = nrow(linpreds))
  }
  
  # Validate parameter matrix column consistency
  param_matrices <- list(location_matrix, scale_matrix, shape_matrix)
  provided_matrices <- param_matrices[!sapply(param_matrices, is.null)]
  if (length(provided_matrices) > 1) {
    ncols <- unique(sapply(provided_matrices, ncol))
    checkmate::assert_true(length(ncols) == 1,
      .var.name = "parameter matrices must have same number of columns")
  }
  
  # Family-specific parameter matrix validation for distributional families
  if (is_distributional_family(family)) {
    required_params <- family$nlp
    provided_matrices <- list(location_matrix, scale_matrix, shape_matrix)
    provided_count <- sum(!sapply(provided_matrices, is.null))
    
    # Check if any parameter matrices are provided (mixed mode validation)
    if (provided_count > 0) {
      # Get expected parameter names dynamically based on family$nlp
      expected_params <- c("location_matrix", "scale_matrix", "shape_matrix")[seq_len(required_params)]
      provided_params <- c("location_matrix", "scale_matrix", "shape_matrix")[!sapply(provided_matrices, is.null)]
      
      # Validate required parameters are provided
      missing_params <- setdiff(expected_params, provided_params)
      if (length(missing_params) > 0) {
        stop(insight::format_error(
          paste0("Family {.field ", family$family, "} requires ",
                 paste(paste0("{.field ", missing_params, "}"), collapse = ", "), ". ",
                 "Provide all required parameter matrices or none.")
        ))
      }
      
      # Validate no extra parameters are provided
      extra_params <- setdiff(provided_params, expected_params)
      if (length(extra_params) > 0) {
        stop(insight::format_error(
          paste0("Family {.field ", family$family, "} does not use ",
                 paste(paste0("{.field ", extra_params, "}"), collapse = ", "), ".")
        ))
      }
    }
  }
  
  rd_fun <- get_family_rd(object)
  
  if (is_distributional_family(family)) {
    # Check if parameter matrices are provided
    param_matrices <- list(location_matrix, scale_matrix, shape_matrix)
    has_param_matrices <- sum(!sapply(param_matrices, is.null)) > 0
    
    if (has_param_matrices) {
      # Reason: parameter matrices contain fitted values, skip lpi extraction and inverse linking
      par_names <- c("location", "scale", "shape")[seq_len(family$nlp)]
      fitted_parameters <- param_matrices[seq_len(family$nlp)]
      names(fitted_parameters) <- par_names
      
      # Validate fitted parameter structure
      checkmate::assert_true(length(fitted_parameters) == family$nlp)
    } else {
      # Reason: extract linear predictors and apply inverse links for standard flow
      parameter_info <- extract_parameter_info_from_lpmat(linpreds, family)
      par_predictions <- split_linear_predictors_by_lpi(linpreds, parameter_info)
      fitted_parameters <- apply_distributional_inverse_links(par_predictions, family)
    }
    
    fitted_matrix <- do.call(cbind, fitted_parameters)
    weights <- rep(1, nrow(linpreds))
    scale_p <- 1
    
    response_pred_vec <- rd_fun(
      mu = fitted_matrix,
      wt = weights,
      scale = scale_p
    )
    
  } else {
    invlink_fun <- get_family_invlink(object)
    
    scale_p <- object[["sig2"]]
    if (is.null(scale_p)) {
      scale_p <- summary(object)[["dispersion"]]
    }
    
    if (!grepl("tweedie", family[["family"]], ignore.case = TRUE)) {
      scale_p <- rep(scale_p, length(linpreds))
    }
    
    weights <- rep(1, length(linpreds))
    expected_pred_vec <- invlink_fun(eta = as.vector(linpreds))
    
    if (grepl("tweedie", family[["family"]], ignore.case = TRUE)) {
      expected_pred_vec <- pmax(expected_pred_vec, 1e-8)
    }
    
    response_pred_vec <- rd_fun(
      mu = expected_pred_vec,
      wt = weights,
      scale = scale_p
    )
  }
  
  return(matrix(response_pred_vec, nrow = nrow(linpreds)))
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
