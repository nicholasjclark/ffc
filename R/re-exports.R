#' @export
#' @importFrom mgcv s
mgcv::s

#' @export
#' @importFrom mgcv te
mgcv::te

#' @export
#' @importFrom mgcv gam
mgcv::gam

#' @export
#' @importFrom mgcv bam
mgcv::bam

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

#' @importFrom generics forecast
#' @export
generics::forecast

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL
