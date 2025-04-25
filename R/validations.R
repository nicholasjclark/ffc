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

#'@noRd
validate_proportional = function(x) {
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
