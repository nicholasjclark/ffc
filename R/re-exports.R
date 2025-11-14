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

#' @export
#' @importFrom mgcv nb
mgcv::nb

#' @export
#' @importFrom mgcv negbin
mgcv::negbin

#' @export
#' @importFrom mgcv betar
mgcv::betar

#' @export
#' @importFrom mgcv cnorm
mgcv::cnorm

#' @export
#' @importFrom mgcv ocat
mgcv::ocat

#' @export
#' @importFrom mgcv scat
mgcv::scat

#' @export
#' @importFrom mgcv ziP
mgcv::ziP

#' @export
#' @importFrom mgcv multinom
mgcv::multinom

#' @export
#' @importFrom mgcv tw
mgcv::tw

#' @export
#' @importFrom mgcv ldTweedie
mgcv::ldTweedie

#' @export
#' @importFrom mgcv rTweedie
mgcv::rTweedie

#' @export
#' @importFrom mgcv twlss
mgcv::twlss

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

#' @importFrom generics forecast
#' @export
generics::forecast

#' @importFrom generics generate
#' @export
generics::generate

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

#' @importFrom utils tail
#' @importFrom tsibble key_vars index_var
NULL
