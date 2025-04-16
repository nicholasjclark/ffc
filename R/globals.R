#' Prevent R CMD Check notes about missing global variables due to
#' dplyr mutates etc...
#' @noRd
utils::globalVariables(
  c(
    ".basis",
    ".time",
    ".estimate",
    ".realisation",
    ".model",
    "functional_ts",
    "ETS"
  )
)
