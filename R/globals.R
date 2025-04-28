#' Prevent R CMD Check notes about missing global variables due to
#' dplyr mutates etc...
#' @noRd
utils::globalVariables(
  c(
    ".basis",
    ".time",
    ".draw",
    ".evaluation",
    ".pred",
    ".rep",
    ".sim",
    ".estimate",
    ".realisation",
    ".model",
    "functional_ts",
    "ETS",
    "ARIMA",
    "d",
    "D",
    "p",
    "P",
    "Q"
  )
)
