#' @noRd
nlist = function(...) {
  m <- match.call()
  dots <- list(...)
  no_names <- is.null(names(dots))
  has_name <- if (no_names) FALSE else nzchar(names(dots))
  if (all(has_name)) {
    return(dots)
  }
  nms <- as.character(m)[-1]
  if (no_names) {
    names(dots) <- nms
  } else {
    names(dots)[!has_name] <- nms[!has_name]
  }
  dots
}

#' @noRd
`c<-` = function(x, value) {
  c(x, value)
}

#' Expand a tbl_ts to include out of sample rows with NAs
#' @noRd
expand_tbl_ts = function(.data, h) {
  min_time <- min(.data$.time)
  max_time <- max(.data$.time) + h

  .expanded <- .data |>
    tidyr::expand(
      tidyr::nesting(.basis, .realisation),
      .time = min_time:max_time
    )
  .data |>
    as.data.frame() |>
    dplyr::right_join(
      .expanded,
      by = dplyr::join_by(.basis, .time, .realisation)
    ) |>
    dplyr::ungroup()
}

#' Expand a tbl_ts to include out of sample index variables
#' @noRd
make_future_data <- function(.data, h = NULL) {
  n <- fabletools::get_frequencies(h, .data, .auto = "smallest")
  if (length(n) > 1) {
    warn("More than one forecast horizon specified, using the smallest.")
    n <- min(n)
  }
  if (is.null(h)) {
    n <- n * 2
  }

  out <- tsibble::new_data(.data, round(n))
  if (tsibble::index_var(.data) == '.time') {
    out <- out |>
      dplyr::mutate(!!attr(.data, 'time_var') := .time)
  }
  return(out)
}

#' Checking functions written by Michell O'hara Wild and the fable team
#' @noRd
check_gaps <- function(x) {
  if (any(tsibble::has_gaps(x)[[".gaps"]])) {
    stop(
      insight::format_error(
        paste0(
          deparse(substitute(x)),
          " contains implicit gaps in time. You should check your data and convert implicit gaps into explicit missing values using `tsibble::fill_gaps()` if required."
        )
      ),
      call. = FALSE
    )
  }
}

#' @noRd
check_regular <- function(x) {
  if (!tsibble::is_regular(x)) {
    stop(
      insight::format_error(
        paste0(
          deparse(substitute(x)),
          " is an irregular time series, which this model does not support. You should consider if your data can be made regular, and use `tsibble::update_tsibble(",
          deparse(substitute(x)),
          ", regular = TRUE)` if appropriate."
        )
      ),
      call. = FALSE
    )
  }
}

#' @noRd
check_ordered <- function(x) {
  if (!tsibble::is_ordered(x)) {
    stop(
      insight::format_error(
        paste0(
          deparse(substitute(x)),
          " is an unordered time series. To use this model, you first must sort the data in time order using `dplyr::arrange(",
          paste(
            c(deparse(substitute(x)), tsibble::key_vars(x)),
            collapse = ", "
          ),
          ", ",
          tsibble::index_var(x),
          ")`"
        )
      ),
      call. = FALSE
    )
  }
}

#' @noRd
all_tsbl_checks <- function(.data) {
  check_gaps(.data)
  check_regular(.data)
  check_ordered(.data)
  if (NROW(.data) == 0) {
    stop(
      insight::format_error(
        "There is no data to model. Please provide a dataset with at least one observation."
      ),
      call. = FALSE
    )
  }
}

#' Extract response variable names from formula
#'
#' Handles both simple responses (y ~ ...) and cbind responses
#' (cbind(successes, failures) ~ ...) robustly
#'
#' @param formula A formula object
#' @param return_all Logical, if TRUE returns all response variable names for
#'   cbind, if FALSE reconstructs cbind() call (default FALSE)
#' @return Character vector of response variable name(s)
#' @noRd
validate_formula_input <- function(formula_or_list, .var.name = "formula") {
  if (is.list(formula_or_list)) {
    checkmate::assert_list(
      formula_or_list,
      types = "formula",
      min.len = 1,
      .var.name = .var.name
    )
  } else {
    checkmate::assert_class(formula_or_list, "formula", .var.name = .var.name)
  }

  return(formula_or_list)
}

#' Get primary formula from single or list input
#'
#' Extracts the primary (first) formula from either a single formula
#' or list of formulae. Used when only one formula is needed for processing.
#'
#' @param formula_or_list Either a single formula or list of formulae
#' @return Single formula object
#' @noRd
get_primary_formula <- function(formula_or_list) {
  validate_formula_input(formula_or_list)

  if (is.list(formula_or_list)) {
    return(formula_or_list[[1]])
  }

  return(formula_or_list)
}

extract_response_vars <- function(formula, return_all = FALSE) {
  # Get primary formula for response extraction
  formula <- get_primary_formula(formula)

  # Get the left-hand side of formula
  lhs <- rlang::f_lhs(formula)

  if (is.null(lhs)) {
    stop(
      insight::format_error("Formula has no response variable"),
      call. = FALSE
    )
  }

  # Handle different types of left-hand sides
  if (typeof(lhs) == "language") {
    # Check if it's cbind specifically
    if (length(lhs) > 1 && as.character(lhs[[1]]) == "cbind") {
      # For cbind, use as.character to get individual variable names
      resp_terms <- as.character(lhs)
    } else {
      # For other transformations (sqrt, log, etc.), use deparse
      resp_terms <- deparse(lhs)
    }
  } else {
    # For simple symbols, use as.character
    resp_terms <- as.character(lhs)
  }

  if (length(resp_terms) == 1L) {
    response <- resp_terms
  } else {
    if (any(grepl('cbind', resp_terms))) {
      resp_terms <- resp_terms[-grepl('cbind', resp_terms)]

      if (return_all) {
        response <- resp_terms
      } else {
        # Reconstruct cbind() call as string
        response <- paste0('cbind(', paste(resp_terms, collapse = ', '), ')')
      }
    } else {
      response <- resp_terms[1]
    }
  }

  return(response)
}

#' Extract underlying variable names from transformed expressions
#'
#' Detects common transformations and extracts the underlying variable names
#' for validation purposes. Optimized with early exit patterns.
#'
#' @param expr_vars Character vector of potentially transformed expressions
#' @return Character vector of underlying variable names
#' @noRd
extract_base_variables <- function(expr_vars) {
  checkmate::assert_character(expr_vars, min.len = 1, any.missing = FALSE)

  # Validate no empty strings
  if (any(trimws(expr_vars) == "")) {
    stop(
      insight::format_error(
        "Empty or whitespace-only expressions not allowed"
      ),
      call. = FALSE
    )
  }

  # Common transformation patterns (ordered by frequency of use)
  patterns <- c(
    "^sqrt\\s*\\(\\s*([^,)]+)\\s*\\)",
    "^log\\s*\\(\\s*([^,)]+)\\s*(?:,.*)?\\)",
    "^log10\\s*\\(\\s*([^,)]+)\\s*\\)",
    "^log2\\s*\\(\\s*([^,)]+)\\s*\\)",
    "^exp\\s*\\(\\s*([^,)]+)\\s*\\)",
    "^abs\\s*\\(\\s*([^,)]+)\\s*\\)",
    "^scale\\s*\\(\\s*([^,)]+)\\s*(?:,.*)?\\)",
    "^sin\\s*\\(\\s*([^,)]+)\\s*\\)",
    "^cos\\s*\\(\\s*([^,)]+)\\s*\\)",
    "^tan\\s*\\(\\s*([^,)]+)\\s*\\)"
  )

  base_vars <- character(length(expr_vars))

  for (i in seq_along(expr_vars)) {
    expr <- trimws(expr_vars[i])
    base_var <- expr # Default to original

    # Check transformation patterns with early exit
    for (pattern in patterns) {
      match_result <- regexec(pattern, expr)
      if (match_result[[1]][1] != -1) {
        captured <- regmatches(expr, match_result)[[1]]
        if (length(captured) >= 2) {
          candidate <- trimws(gsub("^['\"]|['\"]$", "", captured[2]))

          # Validate simple variable name
          if (grepl("^[a-zA-Z][a-zA-Z0-9_\\.]*$", candidate)) {
            base_var <- candidate
            break
          }
        }
      }
    }

    base_vars[i] <- base_var
  }

  return(base_vars)
}


#' Compute functional coefficient predictions using matrix operations
#'
#' @param interpreted_data Data frame containing fts basis evaluations and time
#' @param functional_fc Data frame with functional forecasts (.time, .basis,
#'   .sim, .realisation, .rep columns)
#' @param time_var Character string naming the time variable in interpreted_data
#' @return Data frame with prediction results by draw and observation
#' @noRd
compute_functional_predictions <- function(
  interpreted_data,
  functional_fc,
  time_var
) {
  # Input validation
  checkmate::assert_data_frame(interpreted_data, min.rows = 1)
  checkmate::assert_data_frame(functional_fc, min.rows = 1)
  checkmate::assert_string(time_var)

  # Ensure row identifiers exist for proper tracking
  interpreted_data <- add_row_identifiers(
    interpreted_data,
    ".row_id",
    validate_only = FALSE
  )

  # Validate time variable exists
  if (!time_var %in% names(interpreted_data)) {
    stop(
      insight::format_error(
        paste0("Time variable {.field {time_var}} not found in data")
      ),
      call. = FALSE
    )
  }

  # Validate functional forecast columns
  required_cols <- c(".time", ".basis", ".sim", ".realisation", ".rep")
  missing_cols <- setdiff(required_cols, names(functional_fc))
  if (length(missing_cols) > 0) {
    stop(
      insight::format_error(
        paste0(
          "Missing required columns in functional forecast: ",
          "{.field {paste(missing_cols, collapse = ', ')}}"
        )
      ),
      call. = FALSE
    )
  }

  # Extract and validate fts basis columns (handles both single and multi-parameter models)
  fts_cols <- grep(
    "^(location_|scale_|shape_|param[0-9]+_)?fts_",
    names(interpreted_data),
    value = TRUE
  )
  if (length(fts_cols) == 0) {
    stop(
      insight::format_error(
        "No functional basis columns found. Expected columns starting with 'fts_' or parameter-prefixed 'fts_'"
      ),
      call. = FALSE
    )
  }

  # Validate dimension compatibility
  unique_basis <- unique(functional_fc$.basis)
  if (!setequal(unique_basis, fts_cols)) {
    stop(
      insight::format_error(
        paste0(
          "Mismatch between basis columns in data (",
          length(fts_cols),
          ") and unique basis functions in coefficients (",
          length(unique_basis),
          ")"
        )
      ),
      call. = FALSE
    )
  }

  # Extract basis evaluation matrix
  basis_matrix <- as.matrix(interpreted_data[, fts_cols])
  colnames(basis_matrix) <- fts_cols

  # Convert functional coefficients to wide format
  coeff_wide <- functional_fc |>
    tidyr::pivot_wider(
      id_cols = c(.time, .realisation, .rep),
      names_from = .basis,
      values_from = .sim,
      names_sort = TRUE
    ) |>
    dplyr::arrange(.time, .realisation, .rep)

  # Validate coefficient matrix structure
  coeff_fts_cols <- intersect(fts_cols, names(coeff_wide))
  if (length(coeff_fts_cols) != length(fts_cols)) {
    stop(
      insight::format_error(
        "Not all basis columns found in coefficient wide format"
      ),
      call. = FALSE
    )
  }

  # Compute predictions by row ID instead of positional indexing
  result_list <- list()
  n_obs <- nrow(interpreted_data)

  for (row_idx in seq_len(n_obs)) {
    # Get row ID for tracking (while maintaining .row interface)
    row_id <- interpreted_data$.row_id[row_idx]
    time_point <- interpreted_data[[time_var]][row_idx]

    # Get basis evaluations for this observation
    basis_vals <- basis_matrix[row_idx, fts_cols, drop = FALSE]

    # Get coefficients for this time point
    time_mask <- coeff_wide$.time == time_point
    if (!any(time_mask)) {
      stop(
        insight::format_error(
          "No coefficient data found for time point {.field {time_point}}"
        ),
        call. = FALSE
      )
    }

    time_coeffs <- as.matrix(coeff_wide[time_mask, fts_cols])
    time_meta <- coeff_wide[time_mask, c(".time", ".realisation", ".rep")]

    # Vectorized element-wise multiplication
    pred_matrix <- sweep(time_coeffs, 2, basis_vals, "*")

    # Create result data frame (maintain .row interface, use row_id value)
    result_df <- cbind(
      data.frame(.row = row_id),
      pred_matrix,
      data.frame(.draw = paste0(time_meta$.realisation, "_", time_meta$.rep))
    )

    result_list[[row_idx]] <- result_df
  }

  # Combine results and ensure proper ordering
  fts_fc <- do.call(rbind, result_list) |>
    dplyr::arrange(.row)

  return(fts_fc)
}

#' Add row identifiers to data for tracking through pipeline
#'
#' Adds a unique .row_id column to track original row positions throughout
#' the forecasting pipeline. This ensures results can be mapped back to their
#' original data rows regardless of intermediate sorting operations.
#'
#' @param data A data frame or tibble to add row identifiers to
#' @param id_col Character name of the row ID column (default ".row_id")
#' @param validate_only Logical, if TRUE only validate existing IDs
#' @return The input data with row IDs, or logical if validate_only=TRUE
#' @noRd
add_row_identifiers <- function(
  data,
  id_col = ".row_id",
  validate_only = FALSE
) {
  # Input validation following package standards
  checkmate::assert_data_frame(data, min.rows = 1)
  checkmate::assert_string(
    id_col,
    min.chars = 1,
    pattern = "^[a-zA-Z_\\.][a-zA-Z0-9_\\.]*$"
  )
  checkmate::assert_logical(validate_only, len = 1)

  # Check if row ID column already exists
  has_id_col <- id_col %in% names(data)

  if (has_id_col) {
    # Validate existing IDs are unique and complete
    existing_ids <- data[[id_col]]
    is_valid <- length(unique(existing_ids)) == nrow(data) &&
      !any(is.na(existing_ids))

    if (validate_only) {
      return(is_valid)
    }

    if (is_valid) {
      return(data)
    }

    # Replace with new IDs using proper warning format
    warning(insight::format_warning(
      paste0("Replacing existing ", id_col, " column with new row identifiers")
    ))
  } else if (validate_only) {
    # Column doesn't exist and we're only validating
    return(FALSE)
  }

  # Add sequential row identifiers
  data[[id_col]] <- seq_len(nrow(data))

  return(data)
}

#' Restore data to original row order using row IDs
#'
#' Reorders data back to match original row ordering based on row IDs.
#' This ensures forecast results are returned in the user's expected order.
#'
#' @param data A data frame with row identifiers
#' @param id_col Character name of the row ID column (default ".row_id")
#' @param drop_id Logical whether to remove ID column after ordering
#' @return The data reordered by row IDs
#' @noRd
restore_original_order <- function(data, id_col = ".row_id", drop_id = FALSE) {
  # Input validation
  checkmate::assert_data_frame(data, min.rows = 1)
  checkmate::assert_string(id_col, min.chars = 1)
  checkmate::assert_logical(drop_id, len = 1)

  # Check row ID column exists
  if (!(id_col %in% names(data))) {
    stop(
      insight::format_error(
        paste0(
          "Row ID column {.field ",
          id_col,
          "} not found. Cannot restore original order."
        )
      ),
      call. = FALSE
    )
  }

  # Sort by row ID to restore original order
  data <- data[order(data[[id_col]]), ]

  # Optionally remove the ID column
  if (drop_id) {
    data[[id_col]] <- NULL
  }

  return(data)
}
