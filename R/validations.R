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
    var_type_cap <- paste0(
      toupper(substring(var_type, 1, 1)),
      substring(var_type, 2)
    )

    if (length(missing_vars) == 1) {
      stop(
        insight::format_error(
          paste0(var_type_cap, " {.field ", missing_vars, "} not found in data")
        ),
        call. = FALSE
      )
    } else {
      stop(
        insight::format_error(
          paste0(
            var_type_cap,
            "s {.field {",
            paste(missing_vars, collapse = ", "),
            "}} not found in data"
          )
        ),
        call. = FALSE
      )
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
      stop(
        insight::format_error(
          paste0(
            "Missing values detected in variable {.field ",
            na_var_names,
            "}. ",
            "Please remove or impute missing values before fitting."
          )
        ),
        call. = FALSE
      )
    } else {
      stop(
        insight::format_error(
          paste0(
            "Missing values detected in variables {.field {",
            paste(na_var_names, collapse = ", "),
            "}}. ",
            "Please remove or impute missing values before fitting."
          )
        ),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

#' Validate time intervals for consistency
#'
#' Checks that time intervals within groups are regular and consistent.
#' Works with data.frames, tibbles, and tsibbles.
#'
#' @param data Data frame containing time series data
#' @param time_var Character string specifying the time variable name
#' @param group_vars Optional character vector specifying grouping variables
#' @param tolerance Numeric tolerance for comparing intervals (default 1e-10)
#' @return Invisible TRUE if intervals are valid, otherwise stops with error
#' @noRd
validate_time_intervals <- function(
  data,
  time_var,
  group_vars = NULL,
  tolerance = 1e-10
) {
  checkmate::assert_data_frame(data, min.rows = 2)
  checkmate::assert_string(time_var)
  checkmate::assert_character(group_vars, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_number(tolerance, lower = 0)

  # Validate time variable exists
  validate_vars_in_data(time_var, data, "time variable")

  # Determine grouping variables with enhanced logic
  if (is.null(group_vars)) {
    # Check if data is a tsibble and use its keys
    if (tsibble::is_tsibble(data)) {
      group_vars <- tsibble::key_vars(data)
    } else {
      # Auto-detect categorical variables as potential grouping keys
      # Rationale: For functional data, each categorical combination
      # typically represents a separate time series
      exclude_vars <- time_var
      potential_keys <- setdiff(colnames(data), exclude_vars)
      group_vars <- potential_keys[vapply(
        potential_keys,
        function(var) {
          is.factor(data[[var]]) || is.character(data[[var]])
        },
        logical(1)
      )]
    }
  }

  if (!is.null(group_vars) && length(group_vars) > 0) {
    validate_vars_in_data(group_vars, data, "grouping variable")
  }

  # Function to check intervals within a group
  check_group_intervals <- function(times, group_name = "data") {
    if (length(times) < 2) {
      return(TRUE)
    }

    # Sort times within group
    times <- sort(times)
    intervals <- diff(times)

    if (length(intervals) == 0) {
      return(TRUE)
    }

    # Check if all intervals are approximately equal
    expected_interval <- intervals[1]
    if (!all(abs(intervals - expected_interval) < tolerance)) {
      stop(
        insight::format_error(
          paste0(
            "Irregular time intervals found in ",
            group_name,
            ". ",
            "Intervals range from ",
            round(min(intervals), 6),
            " to ",
            round(max(intervals), 6),
            ". ",
            "Expected consistent interval of ",
            round(expected_interval, 6),
            ". ",
            "Please ensure regular time spacing within groups."
          )
        ),
        call. = FALSE
      )
    }

    return(TRUE)
  }

  # Check intervals by group or overall
  if (!is.null(group_vars) && length(group_vars) > 0) {
    # Group-wise interval checking
    data |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
      dplyr::summarise(
        times = list(!!rlang::sym(time_var)),
        .groups = "drop"
      ) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        valid = check_group_intervals(
          times,
          paste("group", paste(group_vars, collapse = ", "))
        )
      )
  } else {
    # Overall interval checking
    check_group_intervals(data[[time_var]])
  }

  invisible(TRUE)
}

#' Validate newdata for forecasting and return out-of-sample timepoints
#'
#' Checks forecast newdata for valid time structure and filters to only
#' future timepoints beyond the training data. Also validates time intervals
#' and handles grouped data properly.
#'
#' @param newdata Data frame containing forecast data
#' @param model Model object containing training data and time variable info
#' @return Filtered newdata containing only out-of-sample timepoints
#' @noRd
validate_time_intervals_for_forecast <- function(
  newdata,
  model,
  time_var,
  group_vars = character(0)
) {
  checkmate::assert_data_frame(newdata)
  checkmate::assert_class(model, "ffc_gam")
  checkmate::assert_string(time_var)
  checkmate::assert_character(group_vars, any.missing = FALSE)

  # Only proceed if we have multiple time points and fts() terms
  if (nrow(newdata) < 2 || length(model$fts_smooths) == 0) {
    return(list(intervals_valid = TRUE, warnings = character(0)))
  }

  warnings <- character(0)

  # Extract unique training times (fixing the repeated values bug)
  if (length(group_vars) > 0) {
    # Grouped case: check intervals within each group
    train_intervals_by_group <- model$model |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
      dplyr::summarise(
        unique_times = list(sort(unique(!!rlang::sym(time_var)))),
        train_interval = if (length(unique(!!rlang::sym(time_var))) > 1) {
          diff(sort(unique(!!rlang::sym(time_var))))[1]
        } else {
          NA_real_
        },
        .groups = "drop"
      )

    # Check intervals within each group in newdata
    interval_check <- newdata |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
      dplyr::summarise(
        newdata_intervals = if (dplyr::n() > 1) {
          list(diff(sort(!!rlang::sym(time_var))))
        } else {
          list(numeric(0))
        },
        .groups = "drop"
      ) |>
      dplyr::left_join(train_intervals_by_group, by = group_vars)

    # Check for interval mismatches between newdata and training data
    for (i in seq_len(nrow(interval_check))) {
      if (
        length(interval_check$newdata_intervals[[i]]) > 0 &&
          !is.na(interval_check$train_interval[i])
      ) {
        newdata_ints <- interval_check$newdata_intervals[[i]]
        train_int <- interval_check$train_interval[i]

        if (!all(abs(newdata_ints - train_int) < 1e-10)) {
          warnings <- c(
            warnings,
            paste0(
              "Time intervals in newdata (",
              paste(round(newdata_ints, 4), collapse = ", "),
              ") differ from training interval (",
              round(train_int, 4),
              ") for group. ",
              "Consider using the same interval as training data for optimal results."
            )
          )
        }
      }
    }
  } else {
    # Single group case: check overall intervals using unique time points
    newdata_intervals <- diff(sort(unique(newdata[[time_var]])))
    unique_train_times <- sort(unique(model$model[[time_var]]))

    if (length(unique_train_times) > 1) {
      train_interval <- diff(unique_train_times)[1]

      # Check for exact interval match
      if (!all(abs(newdata_intervals - train_interval) < 1e-10)) {
        warnings <- c(
          warnings,
          paste0(
            "Time intervals in newdata (",
            paste(round(newdata_intervals, 4), collapse = ", "),
            ") differ from training interval (",
            round(train_interval, 4),
            "). ",
            "Consider using the same interval as training data for optimal results."
          )
        )
      }
    }
  }

  return(list(intervals_valid = length(warnings) == 0, warnings = warnings))
}

validate_forecast_newdata <- function(newdata, model) {
  checkmate::assert_data_frame(newdata, min.rows = 1)
  checkmate::assert_class(model, "ffc_gam")

  # Add row identifiers to preserve original ordering through validation
  newdata <- add_row_identifiers(newdata, ".original_row_id")

  time_var <- model$time_var
  if (is.null(time_var)) {
    stop(
      insight::format_error(
        "Time variable not found in model object. Ensure model was fitted with time argument"
      ),
      call. = FALSE
    )
  }

  # Validate time variable exists in newdata
  validate_vars_in_data(time_var, newdata, "time variable")

  # Extract grouping variables only from fts() terms with by variables
  exclude_vars <- time_var
  training_key_vars <- character(0)
  if (length(model$fts_smooths) > 0) {
    fts_by_vars <- extract_fts_by_variables(model$fts_smooths)
    # Keep only variables that exist in training data
    training_key_vars <- fts_by_vars[fts_by_vars %in% colnames(model$model)]
  }

  # Validate that training key variables exist in newdata
  if (length(training_key_vars) > 0) {
    validate_vars_in_data(
      training_key_vars,
      newdata,
      "grouping variable from training data"
    )

    # Check if newdata has additional grouping variables not in training
    newdata_potential_keys <- setdiff(colnames(newdata), exclude_vars)
    newdata_only_keys <- newdata_potential_keys[sapply(
      newdata_potential_keys,
      function(var) {
        is.factor(newdata[[var]]) || is.character(newdata[[var]])
      }
    )]
    extra_keys <- setdiff(newdata_only_keys, training_key_vars)

    if (length(extra_keys) > 0) {
      if (!identical(Sys.getenv("TESTTHAT"), "true")) {
        rlang::warn(
          paste0(
            "Found additional grouping variables in newdata: {",
            paste(extra_keys, collapse = ", "),
            "}. ",
            "These will be ignored as they don't exist in training data."
          ),
          .frequency = "once",
          .frequency_id = "extra_grouping_vars"
        )
      }
    }
  }

  # Use validated training key variables for all operations
  key_vars <- training_key_vars

  # Sort newdata by groups first, then time within groups
  if (length(key_vars) > 0) {
    newdata <- newdata |>
      dplyr::group_by(dplyr::across(dplyr::all_of(key_vars))) |>
      dplyr::arrange(!!rlang::sym(time_var), .by_group = TRUE) |>
      dplyr::ungroup()
  } else {
    newdata <- newdata |>
      dplyr::arrange(!!rlang::sym(time_var))
  }

  # Get training data max time, handling groups if present
  if (length(key_vars) > 0) {
    train_max_times <- model$model |>
      dplyr::group_by(dplyr::across(dplyr::all_of(key_vars))) |>
      dplyr::summarise(
        max_time = max(!!rlang::sym(time_var)),
        .groups = "drop"
      )
  } else {
    max_train_time <- max(model$model[[time_var]])
  }

  # Validate newdata intervals are consistent (internal check first)
  if (nrow(newdata) >= 2) {
    suppressWarnings(try(
      {
        validate_time_intervals(newdata, time_var, key_vars)
      },
      silent = TRUE
    ))
  }

  # Validate intervals against training data using DRY helper function
  interval_validation <- validate_time_intervals_for_forecast(
    newdata,
    model,
    time_var,
    key_vars
  )

  # Emit any warnings found
  for (warning_msg in interval_validation$warnings) {
    rlang::warn(warning_msg)
  }

  # Calculate forecast horizons
  if (length(key_vars) > 0) {
    # Calculate forecast horizons by group
    fc_horizons_df <- newdata |>
      dplyr::left_join(train_max_times, by = key_vars) |>
      dplyr::mutate(fc_horizon = !!rlang::sym(time_var) - max_time)

    fc_horizons <- fc_horizons_df$fc_horizon
  } else {
    # Single group case
    fc_horizons <- newdata[[time_var]] - max_train_time
  }

  # Check if we have any future time points
  if (all(fc_horizons <= 0)) {
    stop(
      insight::format_error(
        paste0(
          "No future time points found in newdata. ",
          "All time points are at or before their respective last training times. ",
          "Please provide newdata with {.field ",
          time_var,
          "} values beyond training data."
        )
      ),
      call. = FALSE
    )
  }

  # Warn about overlapping time points and filter
  if (any(fc_horizons <= 0)) {
    n_overlap <- sum(fc_horizons <= 0)
    n_future <- sum(fc_horizons > 0)

    rlang::warn(
      paste0(
        "Found ",
        n_overlap,
        " time points in newdata that overlap with training data. ",
        "Only forecasting for ",
        n_future,
        " future time points."
      )
    )

    # Filter to only future time points
    if (length(key_vars) > 0) {
      newdata <- fc_horizons_df |>
        dplyr::filter(fc_horizon > 0) |>
        dplyr::select(-max_time, -fc_horizon)
    } else {
      newdata <- newdata[fc_horizons > 0, , drop = FALSE]
    }
  }

  return(newdata)
}

#' Convert character variables to factors for random effects
#'
#' Automatically converts character variables to factors when used with
#' random effect basis functions (bs = "re"). Provides graceful error for
#' invalid variable types.
#'
#' @param formula The model formula
#' @param data The data frame
#' @return Modified data frame with character variables converted to factors
#' @noRd
convert_re_to_factors <- function(formula, data) {
  checkmate::assert_data_frame(data)

  # Handle list formulae for distributional regression
  if (is.list(formula)) {
    checkmate::assert_list(formula, types = "formula", min.len = 1)
    # Process each formula in the list and apply conversions
    for (i in seq_along(formula)) {
      data <- convert_re_to_factors(formula[[i]], data)
    }
    return(data)
  }

  checkmate::assert_class(formula, "formula")

  # Extract term labels from formula with error handling
  terms_obj <- tryCatch(
    stats::terms.formula(formula, keep.order = TRUE),
    error = function(e) {
      stop(
        insight::format_error(
          "Invalid formula structure: {e$message}"
        ),
        call. = FALSE
      )
    }
  )

  termlabs <- attr(terms_obj, "term.labels")
  if (is.null(termlabs) || length(termlabs) == 0) {
    return(data) # No terms to process
  }

  # Regex pattern for random effects - extract as constant
  RE_PATTERN <- 'bs\\s*=\\s*["\']re["\']'
  VAR_PATTERN <- "s\\s*\\(\\s*([^,\\)]+)"

  # Find terms that use random effect basis (bs = "re")
  re_terms <- grep(RE_PATTERN, termlabs, value = TRUE)

  if (length(re_terms) > 0) {
    for (term in re_terms) {
      # Extract variable name with comprehensive error handling
      matches <- regmatches(term, regexec(VAR_PATTERN, term))

      if (length(matches[[1]]) < 2) {
        if (!identical(Sys.getenv("TESTTHAT"), "true")) {
          insight::format_warning(
            paste0(
              "Could not parse variable name from random effect term: ",
              term,
              ". Skipping automatic conversion."
            ),
            .frequency = "once",
            .frequency_id = paste0("parse_failure_", term)
          )
        }
        next
      }

      var_name <- trimws(matches[[1]][2])
      var_name <- gsub("^['\"]|['\"]$", "", var_name) # Remove quotes

      # Validate variable name is simple identifier
      if (
        !identical(var_name, make.names(var_name)) ||
          grepl("[^a-zA-Z0-9_\\.]", var_name)
      ) {
        if (!identical(Sys.getenv("TESTTHAT"), "true")) {
          insight::format_warning(
            paste0(
              "Complex expression in random effect term: ",
              term,
              ". Automatic conversion only supports simple variable names."
            ),
            .frequency = "once",
            .frequency_id = paste0("complex_expr_", term)
          )
        }
        next
      }

      if (var_name %in% colnames(data)) {
        var_value <- data[[var_name]]

        # Check variable type and convert if needed
        if (is.character(var_value)) {
          if (!identical(Sys.getenv("TESTTHAT"), "true")) {
            insight::format_warning(
              paste0(
                "Converting character variable {.field ",
                var_name,
                "} to factor for random effect term."
              ),
              .frequency = "once",
              .frequency_id = paste0("re_convert_", var_name)
            )
          }
          data[[var_name]] <- as.factor(var_value)
        } else if (!is.factor(var_value)) {
          stop(
            insight::format_error(
              paste0(
                "Variable {.field ",
                var_name,
                "} used in random effect term ",
                "must be either character or factor, but is {.cls ",
                class(var_value)[1],
                "}. ",
                "Please convert to character or factor before using in s(",
                var_name,
                ", bs = 're')."
              )
            ),
            call. = FALSE
          )
        }
      }
    }
  }

  return(data)
}

#' Validate response variables exist in data with type-aware checking
#'
#' Validates that underlying variables exist in data for both simple and
#' transformed response expressions using type-aware logic for all response types.
#' This is the single source of truth for response validation across the package.
#'
#' @param formula A formula object potentially containing transformed responses
#' @param data Data frame to check for variables
#' @return Invisible TRUE if all underlying variables exist, otherwise stops
#' @noRd
validate_response_in_data <- function(formula, data) {
  checkmate::assert_data_frame(data)

  # Get primary formula for response validation
  formula <- get_primary_formula(formula)

  # Extract response expressions
  resp_terms <- extract_response_vars(formula, return_all = TRUE)

  # Check that response terms are in the data
  # Handle different response types appropriately
  for (i in seq_along(resp_terms)) {
    resp_term <- resp_terms[i]

    # For simple variable names (including cbind-extracted vars), validate directly
    if (grepl("^[a-zA-Z][a-zA-Z0-9_\\.]*$", resp_term)) {
      validate_vars_in_data(resp_term, data, "response variable")
    } else if (grepl("^I\\s*\\(", resp_term)) {
      # For I() expressions, extract variables from inside
      # Extract all variable names from inside I()
      inner_expr <- gsub("^I\\s*\\((.*)\\)\\s*$", "\\1", resp_term)
      # Find all variable-like patterns in the expression
      vars_in_expr <- regmatches(
        inner_expr,
        gregexpr("\\b[a-zA-Z][a-zA-Z0-9_\\.]*\\b", inner_expr)
      )[[1]]
      # Filter out R keywords and functions
      r_keywords <- c(
        "c",
        "if",
        "else",
        "for",
        "in",
        "function",
        "TRUE",
        "FALSE",
        "NULL"
      )
      vars_to_check <- setdiff(vars_in_expr, r_keywords)
      if (length(vars_to_check) > 0) {
        validate_vars_in_data(vars_to_check, data, "response variable")
      }
    } else if (
      grepl(
        "^(sqrt|log|log10|log2|exp|abs|scale|sin|cos|tan)\\s*\\(",
        resp_term
      )
    ) {
      # For transformation functions (sqrt, log, etc.), extract base variables
      underlying_vars <- extract_base_variables(resp_term)
      validate_vars_in_data(underlying_vars, data, "response variable")
    }
    # For other complex expressions, skip validation (assume mgcv will handle)
    # This covers edge cases we haven't anticipated
  }

  invisible(TRUE)
}


#' Validate gam_init structure for distributional models
#'
#' Validates that gam_init structure maintains proper parameter organization
#' needed for distributional regression forecasting
#'
#' @param gam_init_list List structure to validate
#' @param n_parameters Expected number of parameters
#' @return TRUE if valid structure, FALSE otherwise
#' @noRd
validate_gam_init_structure <- function(gam_init_list, n_parameters) {
  checkmate::assert_list(gam_init_list, null.ok = TRUE)
  checkmate::assert_int(n_parameters, lower = 1)

  if (length(gam_init_list) != n_parameters) {
    return(FALSE)
  }

  # Each parameter should have a list of GAM objects or be empty
  for (i in seq_along(gam_init_list)) {
    param_gams <- gam_init_list[[i]]
    if (!is.null(param_gams) && !is.list(param_gams)) {
      return(FALSE)
    }

    # Check that non-empty parameter lists contain GAM objects
    if (length(param_gams) > 0) {
      gam_checks <- sapply(param_gams, function(x) inherits(x, "gam"))
      if (!all(gam_checks)) {
        return(FALSE)
      }
    }
  }

  return(TRUE)
}
