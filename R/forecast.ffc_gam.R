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

# ============================================================================
# COEFFICIENT SAMPLING
# ============================================================================

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
  tryCatch({
    samples <- mgcv::rmvn(n = n_samples, mu = mu, V = V)

    # Verify the sampling worked correctly
    if (any(is.na(samples)) || any(is.infinite(samples))) {
      stop("Generated samples contain invalid values", call. = FALSE)
    }

    return(samples)

  }, error = function(e) {
    stop("Failed to sample coefficients: ", conditionMessage(e), call. = FALSE)
  })
}

#' Validate sample count parameter
#' @param n_samples Number of samples
#' @noRd
validate_sample_count <- function(n_samples) {
  if (!is.numeric(n_samples) ||
      length(n_samples) != 1 ||
      is.na(n_samples) ||
      n_samples <= 0 ||
      n_samples != floor(n_samples)) {
    stop("'n_samples' must be a positive integer", call. = FALSE)
  }
}

#' Validate variance-covariance matrix
#' @param vcov_mat Variance-covariance matrix
#' @noRd
validate_vcov_matrix <- function(vcov_mat) {
  if (!is.matrix(vcov_mat)) {
    stop("Variance-covariance matrix must be a matrix", call. = FALSE)
  }

  if (nrow(vcov_mat) != ncol(vcov_mat)) {
    stop("Variance-covariance matrix must be square", call. = FALSE)
  }

  # Check if matrix is symmetric (within tolerance)
  if (!isSymmetric(vcov_mat, tol = 1e-10)) {
    stop("Variance-covariance matrix must be symmetric", call. = FALSE)
  }

  # Check if matrix is positive definite
  eigenvals <- eigen(vcov_mat, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigenvals <= 1e-12)) {
    stop("Variance-covariance matrix must be positive definite", call. = FALSE)
  }
}

# ============================================================================
# LINEAR PREDICTOR COMPUTATION
# ============================================================================

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

  # Use improved FTS detection
  if (has_fts_terms(object)) {
    # Only use FTS computation if we can verify it's needed and possible

    # Additional safety checks for FTS computation
    if (is.null(object$time_var)) {
      warning("Model has FTS terms but missing time_var, using basic computation")
      return(compute_basic_linear_predictors(object, newdata, beta_samples))
    }

    # FTS computation
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

#' Validate inputs for linear predictor computation
#' @param object Model object
#' @param newdata New data
#' @param beta_samples Coefficient samples
#' @noRd
validate_linpred_inputs <- function(object, newdata, beta_samples) {
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
}

#' Check if model has functional time series terms
#' @param object Model object
#' @return Logical indicating presence of fts terms
#' @noRd
has_fts_terms <- function(object) {

  # Multiple checks for FTS terms

  # Check 1: gam_init component (original check)
  if (!is.null(object$gam_init) && length(object$gam_init) > 0) {
    return(TRUE)
  }

  # Check 2: Look for fts() in the formula
  if (!is.null(object$formula)) {
    formula_text <- deparse(object$formula)
    if (any(grepl("fts\\(", formula_text))) {
      return(TRUE)
    }
  }

  # Check 3: Look for fts() in terms
  if (!is.null(object$terms)) {
    terms_labels <- attr(object$terms, "term.labels")
    if (any(grepl("fts\\(", terms_labels))) {
      return(TRUE)
    }
  }

  # Check 4: Look in smooth terms
  if (!is.null(object$smooth)) {
    smooth_labels <- sapply(object$smooth, function(x) x$label %||% "")
    if (any(grepl("fts", smooth_labels))) {
      return(TRUE)
    }
  }

  # If none of the above, it's likely a simple model
  return(FALSE)
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
  linpred_matrix <- tryCatch({
    predict(object, newdata = newdata, type = "lpmatrix")
  }, error = function(e) {
    stop("Failed to compute linear predictor matrix: ", conditionMessage(e),
         "\nModel class: ", class(object),
         "\nNewdata columns: ", paste(names(newdata), collapse = ", "),
         call. = FALSE)
  })

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
  model_offset <- attr(linpred_matrix, "model.offset")
  if (!is.null(model_offset)) {
    if (length(model_offset) == 1) {
      full_linpreds <- full_linpreds + model_offset
    } else if (length(model_offset) == ncol(full_linpreds)) {
      full_linpreds <- sweep(full_linpreds, 2, model_offset, "+")
    } else {
      warning("Model offset length (", length(model_offset),
              ") does not match number of predictions (", ncol(full_linpreds), ")")
    }
  }

  return(full_linpreds)
}

#' Get linear predictor matrix from model
#' @param object Model object
#' @param newdata New data
#' @return Linear predictor matrix
#' @noRd
get_linear_predictor_matrix <- function(object, newdata) {
  tryCatch({
    predict(object, newdata = newdata, type = "lpmatrix")
  }, error = function(e) {
    stop("Failed to compute linear predictor matrix: ", conditionMessage(e),
         call. = FALSE)
  })
}

#' Validate matrix dimensions for computation
#' @param linpred_matrix Linear predictor matrix
#' @param newdata New data
#' @param beta_samples Coefficient samples
#' @noRd
validate_matrix_dimensions <- function(linpred_matrix, newdata, beta_samples) {
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
}

#' Compute matrix product for linear predictors
#' @param beta_samples Coefficient samples (n_samples x n_coef)
#' @param linpred_matrix Linear predictor matrix (n_pred x n_coef)
#' @return Matrix product (n_samples x n_pred)
#' @noRd
compute_matrix_product <- function(beta_samples, linpred_matrix) {
  # beta_samples is (n_samples x n_coef)
  # linpred_matrix is (n_pred x n_coef)
  # Result should be (n_samples x n_pred)
  tcrossprod(beta_samples, linpred_matrix)
}

#' Add model offset to linear predictors
#' @param linpreds Linear predictors matrix
#' @param linpred_matrix Original linear predictor matrix with offset attribute
#' @return Modified linear predictors matrix
#' @noRd
add_model_offset <- function(linpreds, linpred_matrix) {
  model_offset <- attr(linpred_matrix, "model.offset")

  if (is.null(model_offset)) {
    return(linpreds)
  }

  if (length(model_offset) == 1) {
    # Single offset value - add to all predictions
    return(linpreds + model_offset)
  } else if (length(model_offset) == ncol(linpreds)) {
    # Vector of offsets - add to each column
    return(sweep(linpreds, 2, model_offset, "+"))
  } else {
    warning("Model offset length (", length(model_offset),
            ") does not match number of predictions (", ncol(linpreds), ")")
    return(linpreds)
  }
}

# ============================================================================
# FUNCTIONAL TIME SERIES LINEAR PREDICTORS
# ============================================================================

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

#' Validate requirements for FTS forecasting
#' @param object Model object
#' @noRd
validate_fts_requirements <- function(object) {
  if (is.null(object$time_var)) {
    stop("Model object missing time_var for FTS forecasting", call. = FALSE)
  }
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
  linpreds <- compute_matrix_product(beta_samples, linpred_matrix)

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
add_fts_forecasts <- function(object, newdata, base_linpreds, model, stationary, ...) {

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

#' Check availability of required packages for FTS forecasting
#' @noRd
check_fts_packages <- function() {
  required_packages <- c("fabletools", "tsibble")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop("Packages required for FTS forecasting not available: ",
         paste(missing_packages, collapse = ", "), call. = FALSE)
  }
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

#' Preserve FTS attributes
#' @param functional_coefs New functional coefficients
#' @param intermed_coefs Original intermediate coefficients
#' @return Functional coefficients with preserved attributes
#' @noRd
preserve_fts_attributes <- function(functional_coefs, intermed_coefs) {
  # Add index information if needed
  if (!is.null(attr(intermed_coefs, "index"))) {
    if (!attr(intermed_coefs, "index") %in% names(functional_coefs)) {
      functional_coefs <- functional_coefs %>%
        dplyr::left_join(
          intermed_coefs %>%
            dplyr::select(.time, !!attr(intermed_coefs, "index")) %>%
            dplyr::distinct(),
          by = dplyr::join_by(.time)
        )
    }
  }

  # Set class and attributes
  structure(
    functional_coefs,
    class = c("fts_ts", "tbl_df", "tbl", "data.frame"),
    time_var = attr(intermed_coefs, "time_var"),
    index = attr(intermed_coefs, "index"),
    index2 = attr(intermed_coefs, "index2"),
    interval = attr(intermed_coefs, "interval"),
    summarized = attr(intermed_coefs, "summarized")
  )
}

#' Forecast functional coefficients
#' @param functional_coefs Functional coefficients
#' @param model Forecasting model
#' @param stationary Stationarity constraint
#' @param max_horizon Maximum horizon
#' @param time_var Time variable
#' @param ... Additional arguments
#' @return Forecasted coefficients
#' @noRd
forecast_functional_coefficients <- function(functional_coefs, model, stationary,
                                             max_horizon, time_var, ...) {

  if (model %in% c('ARDF', 'GPDF', 'VARDF') &&
      any(grepl('_mean', unique(functional_coefs$.basis)))) {

    # Forecast non-mean and mean separately
    forecast_mixed_coefficients(functional_coefs, max_horizon, time_var, model, stationary, ...)

  } else {
    # Standard forecasting
    suppressWarnings(
      forecast(
        object = functional_coefs,
        h = max_horizon,
        times = 200,
        model = model,
        stationary = stationary,
        ...
      )
    )
  }
}

#' Forecast mixed coefficient types (mean and non-mean)
#' @param functional_coefs Functional coefficients
#' @param max_horizon Maximum horizon
#' @param time_var Time variable
#' @param model Forecasting model
#' @param stationary Stationarity constraint
#' @param ... Additional arguments
#' @return Mixed forecasts
#' @noRd
forecast_mixed_coefficients <- function(functional_coefs, max_horizon, time_var,
                                        model, stationary, ...) {

  # Forecast non-mean basis functions
  fc_others <- suppressWarnings(
    forecast(
      object = functional_coefs %>% dplyr::filter(!grepl('_mean', .basis)),
      h = max_horizon,
      times = 200,
      model = model,
      stationary = stationary,
      ...
    )
  )

  if ('.time' %in% colnames(fc_others)) {
    fc_others <- fc_others %>% dplyr::select(-.time)
  }

  # Forecast mean basis with ARIMA
  fc_mean <- suppressWarnings(
    forecast(
      object = functional_coefs %>% dplyr::filter(grepl('_mean', .basis)),
      h = max_horizon,
      times = 200,
      model = 'ARIMA',
      stationary = FALSE,
      ...
    )
  ) %>%
    tsibble::as_tibble() %>%
    dplyr::select(dplyr::any_of(names(fc_others)))

  # Add time variable if needed
  if (!time_var %in% colnames(fc_mean) &
      is.null(attr(functional_coefs, "index"))) {
    fc_mean <- fc_mean %>% dplyr::mutate({{time_var}} := .time)
  }

  # Combine forecasts
  fc_others %>%
    dplyr::mutate(.rep = as.character(.rep)) %>%
    dplyr::bind_rows(fc_mean)
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

# ============================================================================
# PREDICTION GENERATION
# ============================================================================

#' Generate predictions from linear predictors
#'
#' @param object Model object
#' @param linpreds Matrix of linear predictors
#' @param type Type of prediction ("link", "expected", "response")
#' @return Matrix of predictions
#' @noRd
generate_predictions <- function(object, linpreds, type) {

  switch(type,
         "link" = linpreds,
         "expected" = posterior_epred(object, linpreds),
         "response" = posterior_predict(object, linpreds),
         stop("Invalid prediction type: ", type, call. = FALSE)
  )
}

# ============================================================================
# OUTPUT FORMATTING
# ============================================================================

#' Format forecast output according to user preferences
#'
#' @param predictions Matrix of predictions
#' @param summary Whether to return summary statistics
#' @param robust Whether to use robust statistics
#' @param probs Quantile probabilities
#' @return Formatted forecast output
#' @noRd
format_forecast_output <- function(predictions,
                                   summary = TRUE,
                                   robust = FALSE,
                                   probs = c(0.025, 0.1, 0.9, 0.975)) {

  if (!summary) {
    return(create_distribution_object(predictions))
  }

  # Compute summary statistics
  stats <- compute_summary_statistics(predictions, robust)
  quantiles <- compute_quantiles(predictions, probs)

  # Create output data frame
  create_summary_dataframe(stats, quantiles, probs)
}

#' Create distribution object from predictions
#' @param predictions Matrix of predictions
#' @return Distribution object
#' @noRd
create_distribution_object <- function(predictions) {
  distributional::dist_sample(
    lapply(seq_len(ncol(predictions)), function(i) predictions[, i])
  )
}

#' Compute summary statistics
#' @param predictions Matrix of predictions
#' @param robust Whether to use robust statistics
#' @return List with estimates and errors
#' @noRd
compute_summary_statistics <- function(predictions, robust) {
  if (robust) {
    estimates <- apply(predictions, 2, median, na.rm = TRUE)
    errors <- apply(abs(predictions - estimates), 2, median, na.rm = TRUE)
  } else {
    estimates <- apply(predictions, 2, mean, na.rm = TRUE)
    errors <- apply(predictions, 2, sd, na.rm = TRUE)
  }

  list(estimates = estimates, errors = errors)
}

#' Compute quantiles
#' @param predictions Matrix of predictions
#' @param probs Quantile probabilities
#' @return Matrix of quantiles
#' @noRd
compute_quantiles <- function(predictions, probs) {
  apply(predictions, 2, quantile, probs = probs, na.rm = TRUE)
}

#' Create summary data frame
#' @param stats Summary statistics
#' @param quantiles Quantile matrix
#' @param probs Quantile probabilities
#' @return Formatted data frame
#' @noRd
create_summary_dataframe <- function(stats, quantiles, probs) {
  out <- data.frame(
    .estimate = stats$estimates,
    .error = stats$errors,
    t(quantiles)
  )

  # Set column names
  colnames(out) <- c(
    ".estimate",
    ".error",
    paste0(".q", formatC(100 * probs, format = "f", digits = 1))
  )

  # Convert to tibble
  structure(out, class = c("tbl_df", "tbl", "data.frame"))
}

# ============================================================================
# VARIANCE INFLATION FUNCTIONS (KEPT SEPARATE FOR SPECIFIC USE CASES)
# ============================================================================

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
