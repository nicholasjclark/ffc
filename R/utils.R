#' @noRd
nlist = function(...) {
  m <- match.call()
  dots <- list(...)
  no_names <- is.null(names(dots))
  has_name <- if (no_names) FALSE else nzchar(names(dots))
  if (all(has_name)) return(dots)
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
expand_tbl_ts = function(.data, h){
  min_time <- min(.data$.time)
  max_time <- max(.data$.time) + h

  .expanded <- .data |>
    tidyr::expand(tidyr::nesting(.basis, .realisation),
                  .time = min_time:max_time)
  .data |>
    as.data.frame() |>
    dplyr::right_join(.expanded,
                      by = dplyr::join_by(.basis, .time, .realisation)) |>
    dplyr::ungroup()
}

#' Expand a tbl_ts to include out of sample index variables
#' @noRd
make_future_data <- function(.data, h = NULL){
  n <- fabletools::get_frequencies(h, .data, .auto = "smallest")
  if(length(n) > 1){
    warn("More than one forecast horizon specified, using the smallest.")
    n <- min(n)
  }
  if(is.null(h)) n <- n*2

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
extract_response_vars <- function(formula, return_all = FALSE) {
  # Get the left-hand side of formula
  lhs <- rlang::f_lhs(formula)
  
  if (is.null(lhs)) {
    stop(insight::format_error("Formula has no response variable"), 
         call. = FALSE)
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
  checkmate::assert_character(expr_vars, min.len = 1, 
                             any.missing = FALSE)
  
  # Validate no empty strings
  if (any(trimws(expr_vars) == "")) {
    stop(insight::format_error(
      "Empty or whitespace-only expressions not allowed"
    ), call. = FALSE)
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
    base_var <- expr  # Default to original
    
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

#' Optimized FC matrix reshaping
#'
#' Reshapes forecast linear predictors into matrix format using vectorized operations
#'
#' @param fc_linpreds Vector of forecast linear predictors
#' @param draw_ids Vector of draw identifiers corresponding to fc_linpreds
#' @param unique_draws Vector of unique draw identifiers in desired order
#' @return Matrix with rows = draws, columns = prediction points
#' @noRd
optimized_fc_matrix_reshape <- function(fc_linpreds, draw_ids, unique_draws) {
  # Input validation following package standards
  checkmate::assert_numeric(fc_linpreds)
  checkmate::assert_vector(draw_ids)
  checkmate::assert_vector(unique_draws)
  checkmate::assert_count(length(unique_draws), .var.name = "number of draws")
  
  n_draws <- length(unique_draws)
  n_cols <- length(fc_linpreds) / n_draws
  
  # Enhanced dimension validation
  if (length(fc_linpreds) %% n_draws != 0) {
    stop(insight::format_error(
      paste0("Dimension mismatch: fc_linpreds length (", 
             length(fc_linpreds), ") not divisible by draws (", 
             n_draws, ")")
    ))
  }
  
  # Vectorized reshaping approach
  draw_factor <- factor(draw_ids, levels = unique_draws)
  draw_indices <- split(seq_along(fc_linpreds), draw_factor)
  
  fc_matrix <- matrix(
    fc_linpreds[unlist(draw_indices, use.names = FALSE)],
    nrow = n_draws,
    ncol = n_cols,
    byrow = TRUE
  )
  
  return(fc_matrix)
}

#' Compute functional coefficient predictions using matrix operations
#'
#' @param interpreted_data Data frame containing fts basis evaluations and time
#' @param functional_fc Data frame with functional forecasts (.time, .basis, 
#'   .sim, .realisation, .rep columns)
#' @param time_var Character string naming the time variable in interpreted_data
#' @return Data frame with prediction results by draw and observation
#' @noRd
compute_functional_predictions <- function(interpreted_data, functional_fc, 
                                          time_var) {
  # Input validation
  checkmate::assert_data_frame(interpreted_data, min.rows = 1)
  checkmate::assert_data_frame(functional_fc, min.rows = 1)
  checkmate::assert_string(time_var)
  
  # Validate time variable exists
  if (!time_var %in% names(interpreted_data)) {
    stop(insight::format_error(
      paste0("Time variable {.field {time_var}} not found in data")
    ), call. = FALSE)
  }
  
  # Validate functional forecast columns
  required_cols <- c(".time", ".basis", ".sim", ".realisation", ".rep")
  missing_cols <- setdiff(required_cols, names(functional_fc))
  if (length(missing_cols) > 0) {
    stop(insight::format_error(
      paste0("Missing required columns in functional forecast: ",
             "{.field {paste(missing_cols, collapse = ', ')}}")
    ), call. = FALSE)
  }
  
  # Extract and validate fts basis columns
  fts_cols <- grep("^fts_", names(interpreted_data), value = TRUE)
  if (length(fts_cols) == 0) {
    stop(insight::format_error(
      "No functional basis columns found. Expected columns starting with 'fts_'"
    ), call. = FALSE)
  }
  
  # Validate dimension compatibility
  unique_basis <- unique(functional_fc$.basis)
  if (!setequal(unique_basis, fts_cols)) {
    stop(insight::format_error(
      paste0("Mismatch between basis columns in data (", length(fts_cols), 
             ") and unique basis functions in coefficients (", 
             length(unique_basis), ")")
    ), call. = FALSE)
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
    stop(insight::format_error(
      "Not all basis columns found in coefficient wide format"
    ), call. = FALSE)
  }
  
  # Compute predictions by time point using matrix operations
  result_list <- list()
  n_obs <- nrow(interpreted_data)
  
  for (time_idx in seq_len(n_obs)) {
    time_point <- interpreted_data[[time_var]][time_idx]
    
    # Get basis evaluations for this observation
    basis_vals <- basis_matrix[time_idx, fts_cols, drop = FALSE]
    
    # Get coefficients for this time point 
    time_mask <- coeff_wide$.time == time_point
    if (!any(time_mask)) {
      stop(insight::format_error(
        paste0("No coefficient data found for time point {time_point}")
      ), call. = FALSE)
    }
    
    time_coeffs <- as.matrix(coeff_wide[time_mask, fts_cols])
    time_meta <- coeff_wide[time_mask, c(".time", ".realisation", ".rep")]
    
    # Vectorized element-wise multiplication
    pred_matrix <- sweep(time_coeffs, 2, basis_vals, "*")
    
    # Create result data frame
    result_df <- cbind(
      data.frame(.row = time_idx),
      pred_matrix,
      data.frame(.draw = paste0(time_meta$.realisation, "_", time_meta$.rep))
    )
    
    result_list[[time_idx]] <- result_df
  }
  
  # Combine results and ensure proper ordering
  fts_fc <- do.call(rbind, result_list) |>
    dplyr::arrange(.row)
  
  return(fts_fc)
}

