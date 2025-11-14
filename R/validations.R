#' Validate that variables exist in data
#'
#' Checks if specified variables exist in a data frame and provides
#' helpful error messages if they don't
#'
#' @param vars Character vector of variable names to check
#' @param data Data frame to check for variables
#' @param var_type Optional string describing the type of variables 
#'   (e.g., "response", "predictor") for better error messages
#' @return Invisible TRUE if all variables exist, otherwise stops with error
#' @noRd
validate_vars_in_data <- function(vars, data, var_type = "variable") {
  checkmate::assert_character(vars, min.len = 1)
  checkmate::assert_data_frame(data)
  checkmate::assert_string(var_type)
  
  missing_vars <- vars[!vars %in% names(data)]
  
  if (length(missing_vars) > 0) {
    # Capitalize first letter of var_type
    var_type_cap <- paste0(toupper(substring(var_type, 1, 1)), 
                          substring(var_type, 2))
    
    if (length(missing_vars) == 1) {
      stop(insight::format_error(
        paste0(var_type_cap, " {.field ", missing_vars, "} not found in data")
      ), call. = FALSE)
    } else {
      stop(insight::format_error(
        paste0(var_type_cap, "s {.field {", paste(missing_vars, collapse = ", "), "}} not found in data")
      ), call. = FALSE)
    }
  }
  
  invisible(TRUE)
}

#' Validate that data contains no missing values
#'
#' Checks if any variables in a data frame contain missing values and provides
#' helpful error messages identifying which variables have NAs
#'
#' @param data Data frame to check for missing values
#' @return Invisible TRUE if no missing values found, otherwise stops with error
#' @noRd
validate_no_missing_values <- function(data) {
  checkmate::assert_data_frame(data)
  
  if (any(is.na(data))) {
    na_vars <- sapply(data, function(x) any(is.na(x)))
    na_var_names <- names(na_vars)[na_vars]
    
    if (length(na_var_names) == 1) {
      stop(insight::format_error(
        paste0("Missing values detected in variable {.field ", na_var_names, "}. ",
               "Please remove or impute missing values before fitting.")
      ), call. = FALSE)
    } else {
      stop(insight::format_error(
        paste0("Missing values detected in variables {.field {", 
               paste(na_var_names, collapse = ", "), "}}. ",
               "Please remove or impute missing values before fitting.")
      ), call. = FALSE)
    }
  }
  
  invisible(TRUE)
}