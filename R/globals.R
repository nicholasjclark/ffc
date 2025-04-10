#' Prevent R CMD Check notes about missing global variables due to
#' dplyr mutates etc...
#' @noRd
utils::globalVariables(c(
  ".basis",
  ".time",
  "functional_ts"
))
