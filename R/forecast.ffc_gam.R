#' Forecasting `ffc_gam` models
#'
#' @importFrom stats median quantile
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
#' @param ... Other arguments to pass to the Stan dynamic factor models
#' @details Computes forecast distributions from fitted `ffc_gam` objects
#' @seealso [ffc_gam()], [fts()], [forecast.fts_ts()]
#' @return Predicted values on the appropriate scale.
#' If `summary == FALSE`, the output is a matrix. If `summary == TRUE`, the output
#' is a tidy `tbl_df / data.frame`
#' @author Nicholas J Clark
#' @export
forecast.ffc_gam <- function(object,
                             newdata,
                             type = "response",
                             model = "ARIMA",
                             stationary = FALSE,
                             summary = TRUE,
                             robust = FALSE,
                             probs = c(0.025, 0.1, 0.9, 0.975),
                             ...) {

  # Validate inputs
  params <- validate_forecast_inputs(
    object = object,
    newdata = newdata,
    type = type,
    model = model,
    stationary = stationary,
    summary = summary,
    robust = robust,
    probs = probs,
    ...
  )

  # Sample coefficients
  beta_samples <- sample_coefficients(object, n_samples = 200)

  # Compute linear predictors
  linpreds <- compute_linear_predictors(
    object = object,
    newdata = newdata,
    beta_samples = beta_samples,
    model = model,
    stationary = stationary,
    ...
  )

  # Generate predictions
  predictions <- generate_predictions(object, linpreds, type)

  # Format output
  format_forecast_output(
    predictions = predictions,
    summary = summary,
    robust = robust,
    probs = probs
  )
}

#' Sample coefficients from model posterior distribution
#'
#' @param object An ffc_gam model object
#' @param n_samples Number of posterior samples to draw
#' @return Matrix of coefficient samples (n_samples x n_coefficients)
#' @noRd
sample_coefficients <- function(object, n_samples = 200) {
  # Input validation
  if (!is.numeric(n_samples) || length(n_samples) != 1 ||
      is.na(n_samples) || n_samples <= 0 || n_samples != floor(n_samples)) {
    stop("'n_samples' must be a positive integer", call. = FALSE)
  }

  # Extract coefficients and variance-covariance matrix
  mu <- coef(object)
  V <- vcov(object)

  # Check for problematic values
  if (any(is.na(mu)) || any(is.infinite(mu))) {
    stop("Model coefficients contain NA or infinite values", call. = FALSE)
  }

  if (any(is.na(V)) || any(is.infinite(V))) {
    stop("Variance-covariance matrix contains NA or infinite values", call. = FALSE)
  }

  # Ensure exact symmetry (fix numerical precision issues)
  V <- (V + t(V)) / 2

  # Check positive definiteness and fix if needed
  eigendecomp <- eigen(V, symmetric = TRUE)
  eigenvals <- eigendecomp$values

  # Replace tiny negative eigenvalues with small positive values
  min_eigenval <- 1e-19
  if (any(eigenvals < min_eigenval)) {
    eigenvals[eigenvals < min_eigenval] <- min_eigenval
    V <- eigendecomp$vectors %*% diag(eigenvals) %*% t(eigendecomp$vectors)
    warning("Adjusted variance-covariance matrix for numerical stability")
  }

  # Sample using the corrected matrix
  samples <- mgcv::rmvn(n = n_samples, mu = mu, V = V)

  # Verify the sampling worked correctly
  if (any(is.na(samples)) || any(is.infinite(samples))) {
    stop("Generated samples contain invalid values", call. = FALSE)
  }

  return(samples)
}

#' Compute linear predictors for forecasting
#'
#' @param object An ffc_gam model object
#' @param newdata Data frame for predictions
#' @param beta_samples Matrix of sampled coefficients
#' @param model Forecasting model type
#' @param stationary Whether to use stationary constraints
#' @param ... Additional arguments
#' @return Matrix of linear predictors (n_samples x n_predictions)
#' @noRd
compute_linear_predictors <- function(object,
                                      newdata,
                                      beta_samples,
                                      model = "ARIMA",
                                      stationary = FALSE,
                                      ...) {

  #  Enhanced input validation
  if (!is.matrix(beta_samples)) {
    stop("'beta_samples' must be a matrix", call. = FALSE)
  }

  if (ncol(beta_samples) != length(coef(object))) {
    stop("'beta_samples' must have ", length(coef(object)),
         " columns to match model coefficients", call. = FALSE)
  }

  if (!is.data.frame(newdata)) {
    stop("'newdata' must be a data.frame", call. = FALSE)
  }

  # Linear predictor computation
  if (has_fts_terms(object)) {

    compute_fts_linear_predictors(
      object = object,
      newdata = newdata,
      beta_samples = beta_samples,
      model = model,
      stationary = stationary,
      ...
    )

  } else {
    # Simple model - use basic computation
    compute_basic_linear_predictors(object, newdata, beta_samples)
  }
}

#' Compute linear predictors for models without fts() terms
#'
#' @param object Model object
#' @param newdata Prediction data
#' @param beta_samples Coefficient samples
#' @return Matrix of linear predictors
#' @noRd
compute_basic_linear_predictors <- function(object, newdata, beta_samples) {

  # Get linear predictor matrix with error handling
  linpred_matrix <- get_linear_predictor_matrix(
    object,
    newdata = newdata
  )

  # Validate dimensions
  if (ncol(linpred_matrix) != ncol(beta_samples)) {
    stop("Dimension mismatch: linear predictor matrix has ", ncol(linpred_matrix),
         " columns but beta_samples has ", ncol(beta_samples), " columns",
         call. = FALSE)
  }

  if (nrow(linpred_matrix) != nrow(newdata)) {
    stop("Dimension mismatch: linear predictor matrix has ", nrow(linpred_matrix),
         " rows but newdata has ", nrow(newdata), " rows",
         call. = FALSE)
  }

  # Compute linear predictors
  # beta_samples: (n_samples x n_coef)
  # linpred_matrix: (n_pred x n_coef)
  # Result: (n_samples x n_pred)
  full_linpreds <- tcrossprod(beta_samples, linpred_matrix)

  # Add model offset if present
  full_linpreds <- add_model_offset(
    full_linpreds,
    linpred_matrix
  )

  return(full_linpreds)
}

#' Compute linear predictors for models with fts() terms
#'
#' @param object Model object
#' @param newdata Prediction data
#' @param beta_samples Coefficient samples
#' @param model Forecasting model
#' @param stationary Stationarity constraint
#' @param ... Additional arguments
#' @return Matrix of linear predictors
#' @noRd
compute_fts_linear_predictors <- function(object,
                                          newdata,
                                          beta_samples,
                                          model,
                                          stationary,
                                          ...) {

  # Validate requirements for FTS forecasting
  validate_fts_requirements(object)

  # Get base linear predictors
  base_linpreds <- get_fts_base_predictors(object, newdata, beta_samples)

  # Add FTS forecasting component
  add_fts_forecasts(
    object = object,
    newdata = newdata,
    base_linpreds = base_linpreds,
    model = model,
    stationary = stationary,
    ...
  )
}

#' Get base linear predictors for FTS models
#' @param object Model object
#' @param newdata New data
#' @param beta_samples Coefficient samples
#' @return Base linear predictors matrix
#' @noRd
get_fts_base_predictors <- function(object, newdata, beta_samples) {
  time_var <- object$time_var

  # Fix time to last training value
  last_train <- max(object$model[[time_var]])
  newdata_shift <- newdata
  newdata_shift[[time_var]] <- last_train

  # Get base predictors
  linpred_matrix <- get_linear_predictor_matrix(object, newdata_shift)
  linpreds <-  tcrossprod(beta_samples, linpred_matrix)

  # Add offset
  add_model_offset(linpreds, linpred_matrix)
}

#' Add FTS forecasts to base predictions
#' @param object Model object
#' @param newdata New data
#' @param base_linpreds Base linear predictors
#' @param model Forecasting model
#' @param stationary Stationarity constraint
#' @param ... Additional arguments
#' @return Enhanced linear predictors
#' @noRd
add_fts_forecasts <- function(
    object,
    newdata,
    base_linpreds,
    model,
    stationary,
    ...) {

  # Check required packages
  check_fts_packages()

  # Get forecast components
  forecast_components <- prepare_fts_forecast_components(
    object = object,
    newdata = newdata,
    model = model,
    stationary = stationary,
    ...
  )

  # Combine with base predictions
  combine_fts_predictions(base_linpreds, forecast_components)
}

#' Prepare FTS forecast components
#' @param object Model object
#' @param newdata New data
#' @param model Forecasting model
#' @param stationary Stationarity constraint
#' @param ... Additional arguments
#' @return Forecast components
#' @noRd
prepare_fts_forecast_components <- function(object, newdata, model, stationary, ...) {
  time_var <- object$time_var

  # Interpret formula and get horizons
  interpreted <- interpret_ffc(
    formula = object$orig_formula,
    data = newdata,
    newdata = newdata,
    gam_init = object$gam_init,
    time_var = time_var
  )

  fc_horizons <- interpreted$data[[time_var]] - max(object$model[[time_var]])
  max_horizon <- max(fc_horizons)

  # Get functional coefficients
  functional_coefs <- prepare_functional_coefficients(object, time_var)

  # Forecast coefficients
  functional_fc <- forecast_functional_coefficients(
    functional_coefs = functional_coefs,
    model = model,
    stationary = stationary,
    max_horizon = max_horizon,
    time_var = time_var,
    ...
  )

  # Compute FTS predictions
  compute_fts_component_predictions(
    interpreted,
    functional_fc,
    time_var,
    attr(functional_coefs, "index"))
}

#' Prepare functional coefficients for forecasting
#' @param object Model object
#' @param time_var Time variable name
#' @return Prepared functional coefficients
#' @noRd
prepare_functional_coefficients <- function(object, time_var) {
  # Extract coefficient time series
  intermed_coefs <- fts_coefs(object, summary = FALSE, times = 200)

  # Calculate mean and shift so last time point is zero
  shift_tail <- function(x) x - tail(x, 1)

  functional_coefs <- intermed_coefs %>%
    dplyr::group_by(.basis, .time) %>%
    dplyr::mutate(.mean = mean(.estimate)) %>%
    dplyr::select(.basis, .time, !!time_var, .mean) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    dplyr::group_by(.basis) %>%
    dplyr::arrange(.time) %>%
    dplyr::mutate(.mean = shift_tail(.mean)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(.realisation = 1, .estimate = .mean) %>%
    tsibble::as_tsibble(key = c(.basis, .realisation), index = .time)

  # Preserve attributes
  preserve_fts_attributes(functional_coefs, intermed_coefs)
}

#' Compute FTS component predictions
#' @param interpreted Interpreted formula data
#' @param functional_fc Functional forecasts
#' @param time_var Time variable name
#' @return FTS component predictions
#' @noRd
compute_fts_component_predictions <- function(
    interpreted,
    functional_fc,
    time_var,
    index_var) {

  # Ensure original time structure is preserved
  if (!is.null(index_var)) {
    # Join original time data using the index attribute
    functional_fc <- functional_fc %>%
      dplyr::left_join(
        interpreted$orig_data %>%
          tibble::as_tibble() %>%
          dplyr::select(!!time_var, !!index_var) %>%
          dplyr::distinct(),
        by = dplyr::join_by(!!index_var)
      ) %>%
      dplyr::rename(.time = !!time_var)

  } else if (time_var %in% colnames(functional_fc)) {
    # Replace .time column with time_var if available
    functional_fc <- functional_fc %>%
      tsibble::as_tibble() %>%
      dplyr::select(-.time) %>%
      dplyr::bind_cols(
        functional_fc %>%
          tsibble::as_tibble() %>%
          dplyr::select(!!time_var) %>%
          dplyr::rename(.time = !!time_var)
      )
  }

  # Join forecasts to basis evaluations
  fts_predictions <- interpreted$data %>%
    dplyr::select(tidyr::matches("^fts_")) %>%
    dplyr::mutate(
      .time = interpreted$data[[time_var]],
      .row = dplyr::row_number()
    ) %>%
    tidyr::pivot_longer(
      cols = tidyr::starts_with("fts_"),
      names_to = ".basis",
      values_to = ".evaluation"
    ) %>%
    dplyr::left_join(
      functional_fc,
      dplyr::join_by(.time, .basis),
      relationship = "many-to-many"
    ) %>%
    dplyr::mutate(.pred = .evaluation * .sim) %>%
    dplyr::select(.basis, .time, .realisation, .rep, .pred, .row) %>%
    tidyr::pivot_wider(names_from = .basis, values_from = .pred) %>%
    dplyr::arrange(.row) %>%
    dplyr::mutate(.draw = paste0(.realisation, "_", .rep)) %>%
    dplyr::select(-.time, -.realisation, -.rep)

  return(fts_predictions)
}

#' Combine FTS predictions with base predictions
#' @param base_linpreds Base linear predictors
#' @param fts_predictions FTS component predictions
#' @return Combined predictions
#' @noRd
combine_fts_predictions <- function(base_linpreds, fts_predictions) {

  # Calculate row-wise predictions for each draw
  fc_linpreds <- fts_predictions %>%
    dplyr::select(-c(.draw, .row)) %>%
    rowSums()

  # Group by draw and create matrix
  unique_draws <- unique(fts_predictions$.draw)
  fts_matrix <- do.call(
    rbind,
    lapply(unique_draws, function(x) {
      fc_linpreds[which(fts_predictions$.draw == x)]
    })
  )

  # Add to base predictions
  base_linpreds + fts_matrix
}

#' Adjust forecast uncertainty using variance inflation
#' @param forecast_df Forecast data frame
#' @param object_sds Object standard deviations
#' @param times Number of time points
#' @param h Forecast horizon
#' @return Adjusted forecast data frame
#' @noRd
adjust_forecast_uncertainty <- function(forecast_df, object_sds, times, h) {
  forecast_df %>%
    dplyr::left_join(object_sds, by = dplyr::join_by(.basis, .realisation)) %>%
    dplyr::mutate(
      .estimate = vctrs::new_vctr(
        purrr::map2(
          .estimate, .sd,
          \(x, y) generate(
            distributional::dist_normal(
              mean = mean(x),
              sd = sqrt(distributional::parameters(x)$sigma^2 + y^2)
            ),
            times = times
          )
        ),
        class = "distribution"
      )
    ) %>%
    tidyr::unnest(cols = .estimate) %>%
    tidyr::unnest(cols = .estimate) %>%
    dplyr::rename(".sim" = ".estimate") %>%
    dplyr::group_by(.basis, .realisation) %>%
    dplyr::mutate(.rep = rep(1:times, h)) %>%
    dplyr::ungroup()
}

#' Simulate forecasts using fable models
#' @param object_tsbl Time series tibble
#' @param model Forecasting model
#' @param h Forecast horizon
#' @param times Number of simulation times
#' @param stationary Whether to use stationary constraints
#' @param object_sds Optional object standard deviations
#' @return Forecast simulations
#' @noRd
fable_forecast <- function(object_tsbl, model, h, times, stationary = FALSE, object_sds = NULL) {

  # Construct model call
  model_call <- create_fable_model_call(model, stationary)

  # Generate forecasts
  if (!is.null(object_sds)) {
    create_adjusted_fable_forecast(object_tsbl, model_call, h, times, object_sds, model)
  } else {
    create_standard_fable_forecast(object_tsbl, model_call, h, times, model)
  }
}

#' Create fable model call
#' @param model Model name
#' @param stationary Whether to use stationary constraints
#' @return Model call
#' @noRd
create_fable_model_call <- function(model, stationary) {
  if (stationary) {
    do.call(`::`, list("fable", model))(
      .estimate,
      order_constraint = (p + q + P + Q <= 6) & (d + D == 0)
    )
  } else {
    do.call(`::`, list("fable", model))(.estimate)
  }
}

#' Create adjusted fable forecast
#' @param object_tsbl Time series tibble
#' @param model_call Model call
#' @param h Horizon
#' @param times Number of times
#' @param object_sds Standard deviations
#' @param model Model name
#' @return Adjusted forecast
#' @noRd
create_adjusted_fable_forecast <- function(object_tsbl, model_call, h, times, object_sds, model) {
  forecast_df <- object_tsbl %>%
    fabletools::model(model_call) %>%
    forecast(h = h) %>%
    dplyr::mutate(.model = model)

  adjust_forecast_uncertainty(
    forecast_df = forecast_df,
    object_sds = object_sds,
    times = times,
    h = h
  )
}

#' Create standard fable forecast
#' @param object_tsbl Time series tibble
#' @param model_call Model call
#' @param h Horizon
#' @param times Number of times
#' @param model Model name
#' @return Standard forecast
#' @noRd
create_standard_fable_forecast <- function(object_tsbl, model_call, h, times, model) {
  object_tsbl %>%
    fabletools::model(model_call) %>%
    fabletools::generate(h = h, times = times) %>%
    dplyr::mutate(.model = model)
}
