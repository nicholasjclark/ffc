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

#' @noRd
check_gaps <- function(x) {
  if (any(tsibble::has_gaps(x)[[".gaps"]])) {
    abort(sprintf("%s contains implicit gaps in time. You should check your data and convert implicit gaps into explicit missing values using `tsibble::fill_gaps()` if required.", deparse(substitute(x))))
  }
}

#' @noRd
check_regular <- function(x) {
  if (!tsibble::is_regular(x)) {
    abort(sprintf("%s is an irregular time series, which this model does not support. You should consider if your data can be made regular, and use `tsibble::update_tsibble(%s, regular = TRUE)` if appropriate.", deparse(substitute(x)), deparse(substitute(x))))
  }
}

#' @noRd
check_ordered <- function(x) {
  if (!tsibble::is_ordered(x)) {
    abort(sprintf(
      "%s is an unordered time series. To use this model, you first must sort the data in time order using `dplyr::arrange(%s, %s)`",
      deparse(substitute(x)), paste(c(deparse(substitute(x)), key_vars(x)), collapse = ", "), index_var(x)
    ))
  }
}

#' @noRd
all_tsbl_checks <- function(.data) {
  check_gaps(.data)
  check_regular(.data)
  check_ordered(.data)
  if (NROW(.data) == 0) {
    abort("There is no data to model. Please provide a dataset with at least one observation.")
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

#' Validate requirements for FTS forecasting
#' @param object Model object
#' @noRd
validate_fts_requirements <- function(object) {
  if (is.null(object$time_var)) {
    stop("Model object missing time_var for FTS forecasting", call. = FALSE)
  }
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

#' Validate forecast inputs and standardize parameters
#'
#' @param object An object of class 'ffc_gam'
#' @param newdata Data frame containing prediction variables
#' @param type Character string specifying prediction type
#' @param model Character string specifying forecasting model
#' @param stationary Logical indicating whether to enforce stationarity
#' @param summary Logical indicating whether to return summary statistics
#' @param robust Logical indicating whether to use robust statistics
#' @param probs Numeric vector of quantile probabilities
#' @param ... Additional arguments
#' @return List of validated parameters
#' @noRd
validate_forecast_inputs <- function(object,
                                     newdata,
                                     type = "response",
                                     model = "ARIMA",
                                     stationary = FALSE,
                                     summary = TRUE,
                                     robust = FALSE,
                                     probs = c(0.025, 0.1, 0.9, 0.975),
                                     ...) {

  # Validate object class
  if (!inherits(object, "ffc_gam")) {
    stop("'object' must be of class 'ffc_gam'", call. = FALSE)
  }

  # Validate newdata
  if (missing(newdata)) {
    stop("'newdata' argument is required", call. = FALSE)
  }

  if (!is.data.frame(newdata)) {
    if (is.list(newdata)) {
      # Try to convert list to data.frame
      tryCatch({
        newdata <- as.data.frame(newdata)
      }, error = function(e) {
        stop("Cannot convert 'newdata' to data.frame: ", conditionMessage(e),
             call. = FALSE)
      })
    } else {
      stop("'newdata' must be a data.frame or list", call. = FALSE)
    }
  }

  if (nrow(newdata) == 0) {
    stop("'newdata' cannot be empty", call. = FALSE)
  }

  # Validate time variable if object has one
  if (!is.null(object$time_var)) {
    time_var <- object$time_var
    if (!time_var %in% names(newdata)) {
      stop("Time variable '", time_var, "' not found in newdata",
           call. = FALSE)
    }

    # Check time variable is numeric
    if (!is.numeric(newdata[[time_var]])) {
      stop("Time variable '", time_var, "' must be numeric", call. = FALSE)
    }

    # Check for missing values in time variable
    if (any(is.na(newdata[[time_var]]))) {
      stop("Time variable '", time_var, "' cannot contain missing values",
           call. = FALSE)
    }
  }

  # Validate type argument
  type <- match.arg(
    arg = type,
    choices = c("link", "expected", "response")
  )

  # Validate stationary constraint logic
  if (stationary && !identical(model, "ARIMA")) {
    stop("stationary = TRUE only works with model = 'ARIMA'",
         call. = FALSE)
  }

  # Validate logical arguments
  if (!is.logical(summary) || length(summary) != 1) {
    stop("'summary' must be a single logical value", call. = FALSE)
  }

  if (!is.logical(robust) || length(robust) != 1) {
    stop("'robust' must be a single logical value", call. = FALSE)
  }

  if (!is.logical(stationary) || length(stationary) != 1) {
    stop("'stationary' must be a single logical value", call. = FALSE)
  }

  # Validate probability values
  if (!is.numeric(probs) ||
      any(is.na(probs)) ||
      any(probs < 0) ||
      any(probs > 1)) {
    stop("'probs' must contain numeric values between 0 and 1",
         call. = FALSE)
  }

  # Return standardized parameters
  list(
    object = object,
    newdata = newdata,
    type = type,
    model = model,
    stationary = stationary,
    summary = summary,
    robust = robust,
    probs = sort(unique(probs)),  # Sort and remove duplicates
    dots = list(...)
  )
}

#' @noRd
validate_pos_integer <- function(x) {
  s <- substitute(x)
  x <- base::suppressWarnings(as.numeric(x))
  if (length(x) != 1L || anyNA(x)) {
    stop("Argument '", s, "' must be a single numeric value", call. = FALSE)
  }

  if (sign(x) != 1) {
    stop("Argument '", s, "' must be a positive integer", call. = FALSE)
  } else {
    if (x %% 1 != 0) {
      stop("Argument '", s, "' must be a positive integer", call. = FALSE)
    }
  }
}

#' @noRd
validate_proportional <- function(x) {
  s <- substitute(x)
  x <- base::suppressWarnings(as.numeric(x))
  if (length(x) != 1L || anyNA(x)) {
    stop("Argument '", s, "' must be a single numeric value", call. = FALSE)
  }

  if (x < 0 || x > 1) {
    stop(
      "Argument '",
      s,
      "' must be a proportion ranging from 0 to 1, inclusive",
      call. = FALSE
    )
  }
}
