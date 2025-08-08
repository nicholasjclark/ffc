#' Generate predictions from linear predictors
#'
#' @param object Model object
#' @param linpreds Matrix of linear predictors
#' @param type Type of prediction ("link", "expected", "response")
#' @return Matrix of predictions
#' @noRd
generate_predictions <- function(object, linpreds, type) {

  switch(type,
         "link" = linpreds,
         "expected" = posterior_epred(object, linpreds),
         "response" = posterior_predict(object, linpreds),
         stop("Invalid prediction type: ", type, call. = FALSE)
  )
}

#' Posterior expectations
#' @noRd
posterior_epred <- function(object,
                            linpreds) {
  # invlink function
  invlink_fun <- get_family_invlink(object)

  # Compute expectations
  expected_pred_vec <- invlink_fun(
    eta = as.vector(linpreds)
  )

  # Convert back to matrix
  out <- matrix(expected_pred_vec,
                nrow = NROW(linpreds)
  )
  return(out)
}

#' Posterior predictions
#' @noRd
posterior_predict <- function(object,
                              linpreds) {
  # rd function if available
  rd_fun <- get_family_rd(object)

  # invlink function
  invlink_fun <- get_family_invlink(object)

  # Dispersion parameter
  scale_p <- object[["sig2"]]
  if (is.null(scale)) {
    scale_p <- summary(object)[["dispersion"]]
  }

  if (!grepl("tweedie", object$family[["family"]],
             ignore.case = TRUE
  )) {
    scale_p <- rep(scale_p, length(linpreds))
  }

  # weights
  weights <- rep(1, length(linpreds))

  # Compute expectations
  expected_pred_vec <- invlink_fun(
    eta = as.vector(linpreds)
  )

  # Now compute response predictions
  response_pred_vec <- rd_fun(
    mu = expected_pred_vec,
    wt = weights,
    scale = scale_p
  )

  # Convert back to matrix
  out <- matrix(response_pred_vec,
                nrow = NROW(linpreds)
  )
  return(out)
}

#' Functions to enable random number generation from mgcv
#' gam / bam objects. Code is modified from severa internal functions written
#' by Gavin Simpson for the gratia R package, which in turn were modified from
#' original code written by Simon Wood for the mgcv R package

# simulator for tweedie LSS models
#' @importFrom rlang .data
#' @importFrom stats rpois rgamma
#' @importFrom tibble tibble
#' @noRd
rtw <- function(mu, p, phi) {
  if (any(p <= 1 | p >= 2)) {
    stop("'p' must be in interval (1, 2)")
  }
  if (any(phi <= 0)) {
    stop("scale parameter 'phi' must be positive")
  }
  if (any(mu < 0)) {
    stop("mean 'mu' must be non-negative")
  }
  lambda <- mu^(2 - p) / ((2 - p) * phi)
  shape <- (2 - p) / (p - 1)
  scale <- phi * (p - 1) * mu^(p - 1)
  N <- rpois(length(lambda), lambda)
  gs <- rep(scale, N)
  tab <- tibble(
    y = rgamma(gs * 0 + 1, shape = shape, scale = gs),
    lab = rep(seq_along(N), N)
  )
  out <- numeric(length(N))
  out[which(N != 0)] <- tab %>%
    dplyr::group_by(.data$lab) %>%
    dplyr::summarise(summed = sum(.data$y)) %>%
    dplyr::pull(.data$summed)
  out
}

#' converts from theta to power parameter `p` given `a` and `b`
#' @noRd
theta_2_power <- function(theta, a, b) {
  i <- theta > 0
  exp_theta_pos <- exp(-theta[i])
  exp_theta_neg <- exp(theta[!i])
  theta[i] <- (b + a * exp_theta_pos) / (1 + exp_theta_pos)
  theta[!i] <- (b * exp_theta_neg + a) / (1 + exp_theta_neg)
  theta
}

#' extracts the `a` and `b` parameters of the model search over which the power
#' parameter is searched for
#' @noRd
get_tw_ab <- function(family) {
  if (family[["family"]] != "twlss") {
    stop("'model' wasn't fitted with 'twlss()' family.", call. = FALSE)
  }
  rfun <- family$residuals
  a <- get(".a", envir = environment(rfun))
  b <- get(".b", envir = environment(rfun))
  c(a, b)
}

#' @importFrom mgcv fix.family.rd
#' @noRd
fix_family_rd <- function(family, ...) {
  # try to fix up the family used by mgcv to add the $rd component
  # for random deviate sampling

  # try the obvious thing first and see if mgcv::fix.family.rd() already handles
  # family
  fam <- mgcv::fix.family.rd(family)

  # if `family` contains a NULL rd we move on, if it is non-null return early
  # as it doesn't need fixing
  if (!is.null(fam$rd)) {
    return(fam)
  }

  # handle special cases
  fn <- fam[["family"]]

  # handle multivariate normal
  if (identical(fn, "Multivariate normal")) {
    # note: mgcv::mvn is documented to ignore prior weights
    # if we ever need to handle weights to scale V, see this post on CV
    # https://stats.stackexchange.com/a/162885/1390
    rd_mvn <- function(V) {
      function(mu, wt, scale) { # function needs to take wt and scale
        mgcv::rmvn(
          n = nrow(mu),
          mu = mu,
          V = V
        )
      }
    }
    fam$rd <- rd_mvn(solve(crossprod(fam$data$R)))
  }
  if (identical(fn, "twlss")) {
    # this uses some helpers to find the `a` and `b` used during fitting and
    # also to convert what `predict()` etc returns (theta) to power parameter
    rd_twlss <- function(a, b) {
      function(mu, wt, scale) {
        rtw(
          mu = mu[, 1], # fitted(model) for twlss is on response scale!
          p = theta_2_power(theta = mu[, 2], a, b),
          phi = exp(mu[, 3])
        )
      }
    }
    tw_pars <- get_tw_ab(fam)
    fam$rd <- rd_twlss(a = tw_pars[1], b = tw_pars[2])
  }

  # return modified family
  fam
}

#' @importFrom mgcv fix.family.rd
#' @importFrom stats family
#' @noRd
get_family_rd <- function(object) {
  if (inherits(object, "glm")) {
    fam <- family(object) # extract family
  } else {
    fam <- object[["family"]]
  }
  ## mgcv stores data simulation funs in `rd`
  fam <- fix_family_rd(fam)
  if (is.null(fam[["rd"]])) {
    stop("Don't yet know how to simulate from family <",
         fam[["family"]], ">",
         call. = FALSE
    )
  }
  fam[["rd"]]
}

#' @noRd
get_family_invlink <- function(object) {
  if (inherits(object, "glm")) {
    fam <- family(object) # extract family
  } else {
    fam <- object[["family"]]
  }
  ## mgcv stores data simulation funs in `rd`
  fam <- fix_family_rd(fam)
  if (is.null(fam[["rd"]])) {
    stop("Don't yet know how to simulate from family <",
         fam[["family"]], ">",
         call. = FALSE
    )
  }
  fam[["linkinv"]]
}
