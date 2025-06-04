#' Define functions with dynamic coefficients in \pkg{ffc} formulae
#'
#' Set up smooth terms with time-varying (dynamic) coefficients for use in
#' \pkg{ffc} models
#'
#' @inheritParams mgcv::te
#' @param ... 	a list of variables that are the covariates that this smooth is a function of.
#' Transformations whose form depends on the values of the data are best avoided here:
#' e.g. `fts(log(x), z)` is fine, but `fts(I(x / sd(x)), z)` is not.
#' @param mean_only `Logical` indicating whether to only include a single basis function
#' to be modelled as time-varying. This can be helpful if you wish to include a temporal
#' trend that can then be forecasted ahead using appropriate time series models.
#' Default is `FALSE`
#' @param k the dimension(s) of the bases used to represent the smooth term.
#' If not supplied then set either to `10` (if only a single covariate is supplied) or
#' to `5 ^ d`, where `d` is the number of covariates supplied.
#' If supplied as a single number then this basis dimension is used for each basis.
#' If supplied as an array then the elements are the dimensions of the component (marginal)
#' bases of the tensor product.
#' @param time_bs a two letter `character` string indicating the (penalized)
#' smoothing basis to use for the time-varying basis coefficients (eg `"tp"` for thin
#' plate regression spline, `"cr"` for cubic regression spline).
#' see \code{\link[mgcv]{smooth.terms}} for an over view of what is available. It is
#' generally recommended to use one of the doubly-penalized bases (i.e. `"cs"` or `"ts"`)
#' as this helps to ensure the resulting basis function coefficient time series are
#' estimated on an appropriate scale for later forecasting
#' @param time_k the dimension of the bases to be used in the time-varying coefficient
#' smooths (see **Details** below). Arbitrarily set to `10` by default.
#' @param time_m the order of the penalty for the time-varying coefficient smooths
#'  (e.g. `2` for normal cubic spline penalty with 2nd derivatives).
#'  Only some smooth classes use this. The `"ps"` class can use a `2` item `array`
#'  giving the basis and penalty order separately.
#' @param by a `factor` variable of the same dimension as each covariate, used to create
#'  a replicate of the smooth for each factor level.
#' @rdname fts
#' @details \pkg{ffc} will evaluate the basis from these smooths and add the basis functions
#' to the internal model design matrix. Once these basis function predictors are added to the data,
#' a model will then be estimated that allows their coefficients to change through time using
#' terms such as
#' `s(time, by = bfun_1, id = 1, bs = 'cr') + s(time, by = bfun_2, id = 1, bs = 'cr') + ...`.
#' By linking the smoothing parameters using the `id` argument, the time-varying function will
#' be efficiently regularised.
#' @return A class \code{fts.spec} object defining the smooth function
#' to be evaluated and configured for estimating time-varying basis function coefficients
#' @author Nicholas J Clark
#' @export
fts <- function(
    ...,
    mean_only = FALSE,
    k = NA,
    time_k = 10,
    bs = "cr",
    time_bs = "ts",
    m = NA,
    time_m = 2,
    d = NA,
    by = NA,
    xt = NULL,
    pc = NULL) {
  validate_pos_integer(time_k)

  # Terms to be smoothed without evaluation
  vars <- as.list(substitute(list(...)))[-1]

  # Dimension of smoother
  dim <- length(vars)

  # Label
  term <- deparse(vars[[1]], backtick = TRUE)
  if (dim > 1) {
    for (i in 2:dim) {
      term[i] <- deparse(vars[[i]], backtick = TRUE)
    }
  }

  for (i in 1:dim) {
    term[i] <- attr(terms(reformulate(term[i])), "term.labels")
  }

  label <- paste("fts_", term[1], sep = "")
  if (dim > 1) {
    for (i in 2:dim) {
      label <- paste(label, "_", term[i], sep = "")
    }
  }

  # By variable
  by_var <- deparse(substitute(by), backtick = TRUE)

  # Create the smooth object for evaluating the basis functions;
  # warnings and errors should propagate directly from mgcv
  sfun <- match.call()
  args <- match(
    c("k", "bs", "m", "d", "by", "xt", "pc"),
    names(sfun),
    0L
  )
  sfun <- sfun[c(1:(dim + 1), args)]
  if (dim > 1L) {
    sfun[[1L]] <- quote(te)
  } else {
    sfun[[1L]] <- quote(s)
  }

  # Return object, being sure to preserve the mgcv class
  return(
    list(
      call = deparse(substitute(sfun)),
      term = term,
      by = by_var,
      time_bs = time_bs,
      time_k = time_k,
      time_m = time_m,
      label = label,
      mean_only = mean_only
    )
  )
}
