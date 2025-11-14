#' Convert ffc_gam forecasts to fable object
#'
#' @param x An `ffc_gam` object used to generate forecasts
#' @param newdata A `data.frame` containing the forecast data with the
#' response variable
#' @param forecasts Optional pre-computed forecasts from `forecast.ffc_gam()`.
#' If NULL, forecasts will be generated using the remaining arguments
#' @param response Character string specifying the response variable name.
#' If NULL, automatically detected from the model formula
#' @param model A character string representing the forecasting model to use
#' if generating forecasts. Default is "ETS"
#' @param key_vars Optional character vector specifying grouping variables.
#' If NULL, automatically detected from categorical variables in newdata
#' @param ... Additional arguments passed to `forecast.ffc_gam()` if
#' generating forecasts
#' @importFrom checkmate assert_class assert_data_frame assert_string
#' assert_character assert_choice
#' @importFrom insight check_if_installed format_error format_warning
#' @importFrom tsibble as_tsibble
#' @importFrom rlang sym syms
#' @importFrom fabletools as_fable
#' @details Converts forecasting output from `ffc_gam` xs into a properly
#' formatted `fable` x that can be used with `fabletools` functions like
#' `autoplot()`, `accuracy()`, and forecast combination methods
#' @return A `fable` x containing forecast distributions and point
#' estimates
#' @seealso [forecast.ffc_gam()]
#' @author Nicholas J Clark
#' @examples
#' \donttest{
#' # Basic usage with automatic detection
#' library(fable)
#' library(tsibble)
#' library(dplyr)
#'
#' # Prepare tourism data
#' tourism_melb <- tourism %>%
#'   filter(Region == "Melbourne", Purpose == "Visiting") %>%
#'   mutate(quarter = as.numeric(format(Quarter, "%q")),
#'          time = row_number())
#'
#' # Split data
#' train <- tourism_melb %>% slice_head(n = 75)
#' test <- tourism_melb %>% slice_tail(n = 5)
#'
#' # Fit model
#' mod <- ffc_gam(
#'   Trips ~ fts(time, mean_only = TRUE, time_k = 50, time_m = 1) +
#'           fts(quarter, k = 4, time_k = 15, time_m = 1),
#'   time = "time", data = train, family = tw(), engine = "gam"
#' )
#'
#' # Convert to fable with auto-detection
#' fc_fable <- as_fable(mod, newdata = test, model = "ETS")
#'
#' # Use fabletools ecosystem
#' autoplot(fc_fable, train)  # Forecast plot
#' accuracy(fc_fable, test)   # Accuracy metrics
#' hilo(fc_fable, level = c(80, 95))  # Prediction intervals
#'
#' # Distribution summaries
#' fc_fable %>%
#'   summarise(
#'     mean_forecast = mean(.dist),
#'     q25 = quantile(.dist, 0.25),
#'     q75 = quantile(.dist, 0.75)
#'   )
#'
#' # With pre-computed forecasts
#' forecasts <- forecast(mod, newdata = test, model = "ETS", summary = FALSE)
#' fc_fable2 <- as_fable(mod, newdata = test, forecasts = forecasts)
#'
#' # With custom parameters
#' fc_fable3 <- as_fable(
#'   mod,
#'   newdata = test,
#'   model = "ETS",
#'   response = "Trips",
#'   key_vars = c("Region", "State")
#' )
#' }
#' @export
as_fable.ffc_gam <- function(x, newdata, forecasts = NULL,
                             response = NULL, model = "ETS",
                             key_vars = NULL, ...) {

  # Input validation
  checkmate::assert_class(x, "ffc_gam")
  checkmate::assert_data_frame(newdata, min.rows = 1)
  checkmate::assert_string(response, null.ok = TRUE)
  checkmate::assert_string(model)
  checkmate::assert_character(key_vars, null.ok = TRUE)

  if (!is.null(forecasts)) {
    checkmate::assert(
      checkmate::check_class(forecasts, "distribution"),
      checkmate::check_class(forecasts, "dist_sample"),
      checkmate::check_matrix(forecasts),
      checkmate::check_numeric(forecasts)
    )
  }

  # Auto-detect response variable if not provided
  if (is.null(response)) {
    if (is.null(x$formula)) {
      stop(insight::format_error(
        "Cannot auto-detect response variable. Please specify {.field response}"
      ))
    }

    # Use utility function to extract response variables robustly
    response <- extract_response_vars(x$formula, return_all = FALSE)
  }

  # Validate response variables exist in newdata
  # For cbind() responses, extract individual variable names for validation
  response_vars_to_check <- extract_response_vars(x$formula, return_all = TRUE)
  validate_vars_in_data(response_vars_to_check, newdata, "response variable")

  # Detect time variable from original model
  time_var <- x$time_var
  if (is.null(time_var)) {
    stop(insight::format_error(
      "Time variable not found in model object. Ensure model was fitted with time argument"
    ))
  }

  validate_vars_in_data(time_var, newdata, "time variable")

  # Generate forecasts if not provided
  if (is.null(forecasts)) {
    forecasts <- forecast(x, newdata = newdata,
                         model = model, summary = FALSE, ...)
  }

  # Validate forecast dimensions match newdata
  if (is.matrix(forecasts) && ncol(forecasts) != nrow(newdata)) {
    stop(insight::format_error(
      "Forecast dimensions {ncol(forecasts)} do not match newdata rows {nrow(newdata)}"
    ))
  }

  # Prepare fable data structure
  fable_data <- newdata

  # Create distribution column using distributional package
  # Reason: Handle different forecast formats for fabletools compatibility
  if (inherits(forecasts, "distribution")) {
    # Already in distribution format
    fable_data[[".dist"]] <- forecasts
    fable_data[[".mean"]] <- as.numeric(mean(forecasts))
  } else if (inherits(forecasts, "dist_sample")) {
    fable_data[[".dist"]] <- forecasts
    fable_data[[".mean"]] <- sapply(forecasts, mean)
  } else if (is.matrix(forecasts)) {
    forecast_list <- lapply(seq_len(ncol(forecasts)), function(i) {
      forecasts[, i]
    })
    fable_data[[".dist"]] <- distributional::dist_sample(forecast_list)
    fable_data[[".mean"]] <- colMeans(forecasts)
  } else {
    fable_data[[".dist"]] <- distributional::dist_degenerate(
      as.numeric(forecasts)
    )
    fable_data[[".mean"]] <- as.numeric(forecasts)
  }

  # Validate time index for tsibble compatibility
  time_classes <- c("Date", "POSIXct", "yearquarter", "yearmonth",
                    "numeric", "integer")

  if (!any(sapply(time_classes, function(cls) {
    inherits(fable_data[[time_var]], cls)
  }))) {
    # Reason: Convert time to numeric if not compatible with tsibble
    if (is.character(fable_data[[time_var]])) {
      insight::format_warning(
        "Converting character time variable to numeric for tsibble compatibility"
      )
      fable_data[[time_var]] <- as.numeric(
        as.factor(fable_data[[time_var]])
      )
    } else {
      stop(insight::format_error(
        "Time variable {.field {time_var}} must be Date, POSIXct, yearquarter, yearmonth, numeric, or integer"
      ))
    }
  }

  # Identify key variables for tsibble
  if (is.null(key_vars)) {
    exclude_vars <- c(time_var, response, ".dist", ".mean")
    potential_keys <- setdiff(colnames(fable_data), exclude_vars)

    key_vars <- potential_keys[sapply(potential_keys, function(var) {
      is.factor(fable_data[[var]]) || is.character(fable_data[[var]])
    })]
  }

  # Convert to tsibble if not already
  if (!inherits(fable_data, "tbl_ts")) {
    # Check if data already has a time index column (like Quarter)
    # Reason: Handle cases where original data has different time index
    yearquarter_cols <- sapply(fable_data, function(x) {
      inherits(x, c("yearquarter", "yearmonth"))
    })

    if (any(yearquarter_cols)) {
      # Use existing time index if present
      index_col <- names(fable_data)[yearquarter_cols][1]
      if (length(key_vars) > 0) {
        fable_data <- tsibble::as_tsibble(
          fable_data,
          index = !!rlang::sym(index_col),
          key = tidyr::all_of(key_vars)
        )
      } else {
        fable_data <- tsibble::as_tsibble(
          fable_data,
          index = !!rlang::sym(index_col)
        )
      }
    } else {
      # Use time_var from model
      if (length(key_vars) > 0) {
        fable_data <- tsibble::as_tsibble(
          fable_data,
          index = !!rlang::sym(time_var),
          key = tidyr::all_of(key_vars)
        )
      } else {
        fable_data <- tsibble::as_tsibble(
          fable_data,
          index = !!rlang::sym(time_var)
        )
      }
    }
  }

  # Build fable x using public tsibble API
  fc_fable <- tsibble::new_tsibble(
    fable_data,
    response = response,
    dist = ".dist",
    model_cn = ".model",
    class = "fbl_ts"
  )

  # Set dimnames for distribution if missing
  if (is.null(dimnames(fc_fable[[".dist"]]))) {
    dimnames(fc_fable[[".dist"]]) <- list(response)
  }
  if (!identical(response, dimnames(fc_fable[[".dist"]]))) {
    dimnames(fc_fable[[".dist"]]) <- list(response)
  }

  # Add model information
  fc_fable[[".model"]] <- paste0("FFC_", model)

  return(fc_fable)
}
