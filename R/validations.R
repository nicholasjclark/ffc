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
