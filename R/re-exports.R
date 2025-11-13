#' Smooth terms for use in GAMs
#'
#' Constructs smooth terms for inclusion in GAM formulae. Use within 
#' [ffc_gam()] for standard smooths or [fts()] for time-varying smooths.
#'
#' @name s
#' @seealso \code{\link[mgcv]{s}}
#' @examples
#' # Standard smooth in ffc_gam
#' # ffc_gam(y ~ s(x), data = data, time = "time")
#' @export
#' @importFrom mgcv s
mgcv::s

#' Tensor product smooth terms
#'
#' Constructs tensor product smooth terms for multidimensional smoothing.
#' Use within [ffc_gam()] formulae for interactions between variables.
#'
#' @name te
#' @seealso \code{\link[mgcv]{te}}
#' @examples
#' # Tensor product smooth in ffc_gam
#' # ffc_gam(y ~ te(x1, x2), data = data, time = "time")
#' @export
#' @importFrom mgcv te
mgcv::te

#' Generalized Additive Model fitting
#'
#' Fits GAMs using mgcv. In ffc, use [ffc_gam()] instead which extends
#' gam functionality with time-varying coefficients.
#'
#' @name gam
#' @seealso \code{\link[mgcv]{gam}}, [ffc_gam()]
#' @examples
#' # Use ffc_gam() for time-varying functionality
#' # mod <- ffc_gam(y ~ fts(x), data = data, time = "time")
#' @keywords internal
#' @export
#' @importFrom mgcv gam
mgcv::gam

#' Large dataset GAM fitting
#'
#' Fits GAMs for large datasets using mgcv::bam. In ffc, specify 
#' `engine = "bam"` in [ffc_gam()] for large datasets.
#'
#' @name bam
#' @seealso \code{\link[mgcv]{bam}}, [ffc_gam()]
#' @examples
#' # Use bam engine for large datasets
#' # mod <- ffc_gam(y ~ fts(x), data = large_data, time = "time", engine = "bam")
#' @keywords internal
#' @export
#' @importFrom mgcv bam
mgcv::bam

#' Negative binomial family
#'
#' Negative binomial distribution family for count data with overdispersion.
#' Use in [ffc_gam()] `family` argument.
#'
#' @name nb
#' @seealso \code{\link[mgcv]{nb}}
#' @examples
#' # ffc_gam(counts ~ fts(x), data = data, family = nb(), time = "time")
#' @keywords internal
#' @export
#' @importFrom mgcv nb
mgcv::nb

#' Negative binomial family (alternative)
#'
#' Alternative negative binomial family specification.
#' Use in [ffc_gam()] `family` argument.
#'
#' @name negbin
#' @seealso \code{\link[mgcv]{negbin}}
#' @examples
#' # ffc_gam(counts ~ fts(x), data = data, family = negbin(), time = "time")
#' @keywords internal
#' @export
#' @importFrom mgcv negbin
mgcv::negbin

#' Beta regression family
#'
#' Beta distribution family for data on (0,1) interval.
#' Use in [ffc_gam()] `family` argument.
#'
#' @name betar
#' @seealso \code{\link[mgcv]{betar}}
#' @examples
#' # ffc_gam(proportions ~ fts(x), data = data, family = betar(), time = "time")
#' @keywords internal
#' @export
#' @importFrom mgcv betar
mgcv::betar

#' Censored normal family
#'
#' Censored normal distribution family for censored continuous data.
#' Use in [ffc_gam()] `family` argument.
#'
#' @name cnorm
#' @seealso \code{\link[mgcv]{cnorm}}
#' @examples
#' # ffc_gam(censored_y ~ fts(x), data = data, family = cnorm(), time = "time")
#' @keywords internal
#' @export
#' @importFrom mgcv cnorm
mgcv::cnorm

#' Ordered categorical family
#'
#' Ordered categorical distribution family for ordinal response data.
#' Use in [ffc_gam()] `family` argument.
#'
#' @name ocat
#' @seealso \code{\link[mgcv]{ocat}}
#' @examples
#' # ffc_gam(ordered_response ~ fts(x), data = data, family = ocat(), time = "time")
#' @keywords internal
#' @export
#' @importFrom mgcv ocat
mgcv::ocat

#' Scaled Student-t family
#'
#' Scaled Student-t distribution family for robust regression.
#' Use in [ffc_gam()] `family` argument or Stan models.
#'
#' @name scat
#' @seealso \code{\link[mgcv]{scat}}
#' @examples
#' # ffc_gam(y ~ fts(x), data = data, family = scat(), time = "time")
#' @keywords internal
#' @export
#' @importFrom mgcv scat
mgcv::scat

#' Zero-inflated Poisson family
#'
#' Zero-inflated Poisson distribution family for count data with excess zeros.
#' Use in [ffc_gam()] `family` argument.
#'
#' @name ziP
#' @seealso \code{\link[mgcv]{ziP}}
#' @examples
#' # ffc_gam(zero_inflated_counts ~ fts(x), data = data, family = ziP(), time = "time")
#' @keywords internal
#' @export
#' @importFrom mgcv ziP
mgcv::ziP

#' Multinomial family
#'
#' Multinomial distribution family for categorical response data.
#' Use in [ffc_gam()] `family` argument.
#'
#' @name multinom
#' @seealso \code{\link[mgcv]{multinom}}
#' @examples
#' # ffc_gam(categories ~ fts(x), data = data, family = multinom(), time = "time")
#' @keywords internal
#' @export
#' @importFrom mgcv multinom
mgcv::multinom

#' Tweedie family
#'
#' Tweedie distribution family for continuous data with zeros.
#' Use in [ffc_gam()] `family` argument.
#'
#' @name tw
#' @seealso \code{\link[mgcv]{tw}}
#' @examples
#' # ffc_gam(continuous_with_zeros ~ fts(x), data = data, family = tw(), time = "time")
#' @keywords internal
#' @export
#' @importFrom mgcv tw
mgcv::tw

#' Tweedie log density
#'
#' Computes log density for Tweedie distribution.
#' Utility function for Tweedie family.
#'
#' @name ldTweedie
#' @seealso \code{\link[mgcv]{ldTweedie}}, [tw()]
#' @examples
#' # ldTweedie(y = data, mu = fitted_values, p = 1.5, phi = 1)
#' @keywords internal
#' @export
#' @importFrom mgcv ldTweedie
mgcv::ldTweedie

#' Random Tweedie generation
#'
#' Generates random numbers from Tweedie distribution.
#' Utility function for Tweedie family.
#'
#' @name rTweedie
#' @seealso \code{\link[mgcv]{rTweedie}}, [tw()]
#' @examples
#' # rTweedie(n = 100, mu = 1, p = 1.5, phi = 1)
#' @keywords internal
#' @export
#' @importFrom mgcv rTweedie
mgcv::rTweedie

#' Tweedie location-scale-shape family
#'
#' Extended Tweedie family allowing shape and scale parameters to vary.
#' Use in [ffc_gam()] `family` argument for complex Tweedie models.
#'
#' @name twlss
#' @seealso \code{\link[mgcv]{twlss}}, [tw()]
#' @examples
#' # ffc_gam(list(y ~ fts(x), ~ s(z), ~ s(w)), family = twlss(), data = data, time = "time")
#' @keywords internal
#' @export
#' @importFrom mgcv twlss
mgcv::twlss

#' Automatic plotting
#'
#' Automatic plotting method for ffc objects. In ffc, used for plotting
#' time-varying coefficients from [fts_coefs()].
#'
#' @name autoplot
#' @seealso \code{\link[ggplot2]{autoplot}}, [autoplot.fts_ts()]
#' @examples
#' # coefs <- fts_coefs(model, summary = FALSE)
#' # autoplot(coefs)
#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

#' Forecasting methods
#'
#' Generic forecasting function. In ffc, use with [forecast.fts_ts()] and
#' [forecast.ffc_gam()] for functional forecasting.
#'
#' @name forecast
#' @seealso \code{\link[generics]{forecast}}, [forecast.fts_ts()], [forecast.ffc_gam()]
#' @examples
#' # coefs <- fts_coefs(model)
#' # forecast(coefs, h = 5, model = "ARIMA")
#' @importFrom generics forecast
#' @export
generics::forecast

#' Generate method
#'
#' Generic generate function from the generics package.
#' Used internally for forecast simulation.
#'
#' @name generate
#' @seealso \code{\link[generics]{generate}}
#' @examples
#' # Used internally by forecast methods
#' @keywords internal
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

# Import missing functions to fix R CMD check NOTEs
#' @importFrom utils tail
#' @importFrom tsibble key_vars index_var
NULL
