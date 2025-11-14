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
        "{var_type_cap} {.field {missing_vars}} not found in data"
      ), call. = FALSE)
    } else {
      stop(insight::format_error(
        "{var_type_cap}s {.field {missing_vars}} not found in data"
      ), call. = FALSE)
    }
  }
  
  invisible(TRUE)
}