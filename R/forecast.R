#' Forecasting functional basis coefficients
#'
#' @param object An object of class `fts_ts` containing time-varying
#' basis function coefficients extracted from an `ffc_gam` object
#' @param model A `character` string representing a valid univariate model definition
#' from the \pkg{fable} package. Note that the chosen method must have an associated
#' `generate()` method in order to simulate forecast realisations. Valid models
#' currently include: `'ETS'`, `'ARIMA'`, `'AR'`, `'NAIVE'`, and `'NNETAR'`
#' @param h A positive `integer` specifying the length of the forecast
#' horizon
#' @param times A positive `integer` specifying the number of forecast
#' realisation paths to simulate from the fitted forecast `model`
#' @param stationary If `TRUE`, the fitted time series models are
#' constrained to be stationary. Default is `FALSE`. This option
#' only works when `model == 'ARIMA'`
#' @param ... Ignored
#' @details Basis function coefficient time series will be used as input
#' to the specified `model` to train a forecasting model that will then be
#' extrapolated `h` timesteps into the future. A total of `times` forecast realisations
#' will be returned for each basis coefficient
#' @return A `tsibble` object containing the forecast prediction (`.sim`) for each
#' replicate realisation (`.rep`) at each timestep in the forecast horizon `h`
#' @seealso [predict()], [fts()]
#' @author Nicholas J Clark
#' @export
forecast.fts_ts <- function(
    object,
    model = "ARIMA",
    h = 1,
    times = 25,
    stationary = FALSE,
    ...) {
  insight::check_if_installed("fable")

  # Check arguments
  validate_pos_integer(h)
  if (stationary & !identical(model, "ARIMA")) {
    stop("stationary = TRUE only works with model = 'ARIMA'")
  }

  # Convert time-varying coefficients to tsibble
  object_tsbl <- fts_ts_2_tsbl(object)

  # If .sd is included, use the variance inflation method
  # to update forecast uncertainties
  if (exists(".sd", object_tsbl)) {
    # Calculate SDs
    object_sds <- object_tsbl %>%
      as.data.frame() %>%
      dplyr::group_by(.basis, .realisation) %>%
      dplyr::summarise(.sd = mean(.sd))

    if (stationary) {
      out <- object_tsbl %>%
        # Fit the forecasting model(s) to each basis coefficient
        # time series
        fabletools::model(
          do.call(
            what = `::`,
            args = list("fable", model)
          )(.estimate,
            order_constraint = (p + q + P + Q <= 6) & (d + D == 0))
        ) %>%
        # Simulate future paths for each coefficient time series
        forecast(
          h = h
        ) %>%
        # Tidy the names
        dplyr::mutate(.model = model) %>%
        # Add the SDs to the tibble and use them to modify
        # the forecast uncertainty
        dplyr::left_join(object_sds,
          by = dplyr::join_by(.basis, .realisation)
        ) %>%
        dplyr::mutate(.estimate = vctrs::new_vctr(
          purrr::map2(
            .estimate,
            .sd,
            \(x, y) generate(
              distributional::dist_normal(
                mean = mean(x),
                sd = sqrt(distributional::parameters(x)$sigma^2 + y^2)
              ),
              times = times
            )
          ),
          class = "distribution"
        )) %>%
        tidyr::unnest(
          cols = .estimate
        ) %>%
        tidyr::unnest(
          cols = .estimate
        ) %>%
        dplyr::rename(".sim" = ".estimate") %>%
        dplyr::group_by(.basis, .realisation) %>%
        dplyr::mutate(.rep = rep(1:times, h)) %>%
        dplyr::ungroup()
    } else {
      out <- object_tsbl %>%
        # Fit the forecasting model(s) to each basis coefficient
        # time series
        fabletools::model(
          do.call(
            what = `::`,
            args = list("fable", model)
          )(.estimate)
        ) %>%
        # Simulate future paths for each coefficient time series
        forecast(
          h = h
        ) %>%
        # Tidy the names
        dplyr::mutate(.model = model) %>%
        # Add the SDs to the tibble and use them to modify
        # the forecast uncertainty
        dplyr::left_join(object_sds,
          by = dplyr::join_by(.basis, .realisation)
        ) %>%
        dplyr::mutate(.estimate = vctrs::new_vctr(
          purrr::map2(
            .estimate,
            .sd,
            \(x, y) generate(
              distributional::dist_normal(
                mean = mean(x),
                sd = sqrt(distributional::parameters(x)$sigma^2 + y^2)
              ),
              times = times
            )
          ),
          class = "distribution"
        )) %>%
        tidyr::unnest(
          cols = .estimate
        ) %>%
        tidyr::unnest(
          cols = .estimate
        ) %>%
        dplyr::rename(".sim" = ".estimate") %>%
        dplyr::group_by(.basis, .realisation) %>%
        dplyr::mutate(.rep = rep(1:times, h)) %>%
        dplyr::ungroup()
    }
  } else {
    # Else do not modify forecast variances
    if (stationary) {
      out <- object_tsbl %>%
        # Fit the forecasting model(s) to each basis coefficient
        # time series
        fabletools::model(
          do.call(
            what = `::`,
            args = list("fable", model)
          )(.estimate,
            order_constraint = (p + q + P + Q <= 6) & (d + D == 0))
        ) %>%
        # Simulate future paths for each coefficient time series
        fabletools::generate(
          h = h,
          times = times
        ) %>%
        # Tidy the names
        dplyr::mutate(.model = model)
    } else {
      out <- object_tsbl %>%
        fabletools::model(
          do.call(
            what = `::`,
            args = list("fable", model)
          )(.estimate)
        ) %>%
        fabletools::generate(
          h = h,
          times = times
        ) %>%
        dplyr::mutate(.model = model)
    }
  }



  return(out)
}

#' Forecasting `ffc_gam` models
#'
#' @importFrom stats median quantile
#' @inheritParams forecast.fts_ts
#' @param object An object of class `ffc_gam`. See [ffc_gam()]
#' @param newdata `dataframe` or `list` of test data containing the
#' same variables that were included in the original `data` used to fit the
#' model. The covariate information in `newdata`, along with the temporal information
#' indexed by the `time` variable in the original call to `ffc_gam()`, will be used to
#' generate forecasts from the fitted model equations
#' @param type When this has the value `link`, the linear predictor is calculated
#' on the link scale. If `expected` is used, predictions reflect the expectation
#' of the response (the mean) but ignore uncertainty in the observation process.
#' When `response` is used (the default), the predictions take uncertainty in the observation
#' process into account to return predictions on the outcome scale
#' @param quantile_fc Logical. If `TRUE`, quantile forecasts are computed. If
#' `FALSE`, forecasts are computed using a meta-analysis approach that adjusts
#' the forecast variance based on the empirical variance of the basis function
#' coefficient time series. Using `TRUE` will be faster but using `FALSE` may give
#' better forecast uncertainties when the variance of the basis fucntion coefficient
#' time series is not stable over time. Defaults to `FALSE`
#' @param summary Should summary statistics be returned instead of the raw values?
#' Default is `TRUE`
#' @param robust If `FALSE` (the default) the mean is used as the measure of
#' central tendency and the standard deviation as the measure of variability.
#' If `TRUE`, the median and the median absolute deviation (MAD) are applied instead.
#' Only used if `summary` is `TRUE`
#' @param probs The percentiles to be computed by the `quantile()` function.
#' Only used if `summary` is `TRUE`
#' @param ... ignored
#' @details Computes forecast distributions from fitted `ffc_gam` objects
#' @seealso [ffc_gam()], [fts()], [forecast.fts_ts()]
#' @return Predicted values on the appropriate scale.
#' If `summary == FALSE`, the output is a matrix. If `summary == TRUE`, the output
#' is a tidy `tbl_df / data.frame`
#' @author Nicholas J Clark
#' @export
forecast.ffc_gam <- function(
    object,
    newdata,
    type = "response",
    quantile_fc = FALSE,
    model = "ARIMA",
    stationary = FALSE,
    summary = TRUE,
    robust = TRUE,
    probs = c(0.025, 0.1, 0.9, 0.975),
    ...) {
  type <- match.arg(
    arg = type,
    choices = c(
      "link",
      "expected",
      "response"
    )
  )

  validate_proportional(min(probs))
  validate_proportional(max(probs))

  # Extract the full linear predictor matrix
  orig_lpmat <- predict(
    object,
    newdata = newdata,
    type = "lpmatrix"
  )

  # Take full draws of beta coefficients
  orig_betas <- mgcv::rmvn(
    n = 1000,
    mu = coef(object),
    V = vcov(object)
  )

  if (length(object$gam_init) == 0) {
    # No need to modify lpmatrix if there were no
    # fts() terms in the model
    full_linpreds <- matrix(
      as.vector(
        t(apply(
          as.matrix(orig_betas),
          1,
          function(row) {
            orig_lpmat %*% row +
              attr(orig_lpmat, "model.offset")
          }
        ))
      ),
      nrow = NROW(orig_betas)
    )
  } else {
    # Determine horizon (assuming equal time gaps)
    time_var <- object$time_var

    interpreted <- interpret_ffc(
      formula = object$orig_formula,
      data = newdata,
      newdata = newdata,
      gam_init = object$gam_init,
      time_var = object$time_var
    )

    fc_horizons <- interpreted$data[[time_var]] -
      max(object$model[[time_var]])

    max_horizon <- max(fc_horizons)

    # Extract functional basis coefficient time series
    intermed_coefs <- fts_coefs(
      object,
      summary = FALSE,
      times = 1000
    )

    # Calculate man + SD or quantiles
    if (quantile_fc) {
      functional_coefs <- intermed_coefs %>%
        dplyr::group_by(.basis, .time) %>%
        dplyr::summarise(
          tibble::as_tibble_row(
            quantile(
              .estimate,
              probs = seq(0.025, 0.975, length.out = 5),
              type = 8
            )
          )
        ) %>%
        dplyr::ungroup() %>%
        tidyr::pivot_longer(
          cols = tidyr::contains("%"),
          names_to = "quantile",
          values_to = ".estimate"
        ) %>%
        dplyr::mutate(.realisation = as.numeric(gsub("%", "", quantile))) %>%
        tsibble::as_tsibble(
          key = c(.basis, .realisation),
          index = .time
        )
    } else {
      functional_coefs <- intermed_coefs %>%
        dplyr::group_by(.basis, .time) %>%
        dplyr::summarise(
          .mean = mean(.estimate),
          .sd = sd(.estimate)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          .realisation = 1,
          .estimate = .mean
        ) %>%
        tsibble::as_tsibble(
          key = c(.basis, .realisation),
          index = .time
        )
    }


    if (!is.null(attr(intermed_coefs, "index"))) {
      functional_coefs <- functional_coefs %>%
        dplyr::left_join(
          intermed_coefs %>%
            dplyr::select(.time, !!attr(intermed_coefs, "index")) %>%
            dplyr::distinct(),
          by = dplyr::join_by(.time)
        )
    }
    class(functional_coefs) <- c("fts_ts", "tbl_df", "tbl", "data.frame")
    attr(functional_coefs, "time_var") <- attr(intermed_coefs, "time_var")
    attr(functional_coefs, "index") <- attr(intermed_coefs, "index")
    attr(functional_coefs, "index2") <- attr(intermed_coefs, "index2")
    attr(functional_coefs, "interval") <- attr(intermed_coefs, "interval")
    attr(functional_coefs, "summarized") <- attr(intermed_coefs, "summarized")

    # Fit the time series model to the basis coefficients
    # and generate basis coefficient forecasts
    functional_fc <- suppressWarnings(
      forecast(
        object = functional_coefs,
        h = max_horizon,
        times = 1000,
        model = model,
        stationary = stationary,
      )
    )

    if (!is.null(attr(intermed_coefs, "index"))) {
      functional_fc <- functional_fc %>%
        dplyr::left_join(
          interpreted$orig_data %>%
            dplyr::select(!!time_var, !!attr(intermed_coefs, "index")) %>%
            dplyr::distinct(),
          by = dplyr::join_by(!!attr(intermed_coefs, "index"))
        ) %>%
        dplyr::rename(.time = !!time_var)
    }

    # Only need to return forecasts for those times that are
    # in newdata
    functional_fc <- functional_fc %>%
      dplyr::filter(
        .time %in% unique(interpreted$data[[time_var]])
      )

    if (quantile_fc) {
      # Sample draws based on quantile weights
      inds_keep <- functional_fc %>%
        tibble::as_tibble() %>%
        dplyr::select(.basis, .realisation, .rep) %>%
        dplyr::mutate(.weight = norm_quantiles(.realisation)) %>%
        dplyr::group_by(.basis) %>%
        dplyr::sample_n(
          size = 1000,
          replace = TRUE,
          weight = .weight
        )

      functional_fc <- inds_keep %>%
        dplyr::left_join(functional_fc,
          by = dplyr::join_by(.basis, .realisation, .rep),
          relationship = "many-to-many"
        ) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.basis, .time) %>%
        dplyr::mutate(
          .realisation = 1,
          .rep = 1:1000
        )
    }

    # Which coefficients in lpmatrix are associated with fts objects?
    smooth_names <- unlist(
      purrr::map(object$smooth, "label"),
      use.names = FALSE
    )

    fts_names <- grep(
      ":fts_",
      smooth_names,
      fixed = TRUE
    )

    fts_coefs <- unlist(
      purrr::map(
        object$smooth[fts_names],
        \(x) x$first.para:x$last.para
      ),
      use.names = FALSE
    )

    # Drop these columns from the lpmatrix and the betas
    intermed_lpmat <- orig_lpmat[, -fts_coefs, drop = FALSE]
    intermed_betas <- orig_betas[, -fts_coefs, drop = FALSE]

    # Calculate intermediate linear predictors
    intermed_linpreds <- matrix(
      as.vector(
        t(apply(
          as.matrix(intermed_betas),
          1,
          function(row) {
            intermed_lpmat %*% row +
              attr(orig_lpmat, "model.offset")
          }
        ))
      ),
      nrow = NROW(intermed_betas)
    )

    # Join forecasts to the basis function evaluations
    fts_fc <- data.frame(
      interpreted$data[grep("fts_", colnames(interpreted$data))]
    ) %>%
      dplyr::bind_cols(
        data.frame(
          .time = interpreted$data[[time_var]],
          .row = 1:length(interpreted$data[[1]])
        )
      ) %>%
      tidyr::pivot_longer(
        cols = !tidyr::contains(c(
          ".time",
          ".row"
        )),
        names_to = ".basis",
        values_to = ".evaluation"
      ) %>%
      dplyr::left_join(
        functional_fc,
        dplyr::join_by(.time, .basis),
        relationship = "many-to-many"
      ) %>%
      # Calculate prediction
      dplyr::mutate(.pred = .evaluation * .sim) %>%
      # Now take 'draws' of the betas
      dplyr::select(
        .basis, .time, .realisation, .rep, .pred, .row
      ) %>%
      # Pivot back to wide format
      tidyr::pivot_wider(
        names_from = .basis,
        values_from = .pred
      ) %>%
      dplyr::arrange(.row) %>%
      dplyr::mutate(.draw = paste0(.realisation, "_", .rep)) %>%
      dplyr::select(
        -.time,
        -.realisation,
        -.rep
      )

    # Should now have n_draws * n_sims draws for each row of newdata
    # if (!NROW(orig_lpmat) * n_draws * n_sims == NROW(fts_fc)) {
    #   stop("Wrong dimensions in forecast coefs; need to check on this")
    # }

    # Should also have same ncols as number of fts basis functions
    if (
      !NCOL(fts_fc) - 2 ==
        NCOL(interpreted$data[
          grep(
            "fts_",
            colnames(interpreted$data)
          )
        ])
    ) {
      stop("Wrong dimensions in forecast coefs; need to check on this")
    }

    # If dimensions correct, take rowsums for each draw
    fc_linpreds <- fts_fc %>%
      dplyr::select(-c(.draw, .row)) %>%
      rowSums()

    # Add the draw-specific row predictions to the
    # intermediate prediction matrix
    unique_draws <- unique(fts_fc$.draw)
    full_linpreds <- do.call(
      rbind,
      lapply(
        unique_draws,
        function(x) {
          fc_linpreds[which(fts_fc$.draw == x)]
        }
      )
    ) +
      intermed_linpreds
  }

  # Now can proceed to send full_linpreds to the relevant
  # invlink and rng functions for outcome-level predictions
  if (type == "link") preds <- full_linpreds

  if (type == "expected") {
    preds <- posterior_epred(
      object,
      full_linpreds
    )
  }

  if (type == "response") {
    preds <- posterior_predict(
      object,
      full_linpreds
    )
  }

  # Summarise if necessary
  if (summary) {
    Qs <- apply(preds, 2, quantile, probs = probs, na.rm = TRUE)

    if (robust) {
      estimates <- apply(preds, 2, median, na.rm = TRUE)
      errors <- apply(abs(preds - estimates), 2, median, na.rm = TRUE)
    } else {
      estimates <- apply(preds, 2, mean, na.rm = TRUE)
      errors <- apply(preds, 2, sd, na.rm = TRUE)
    }

    out <- data.frame(
      cbind(estimates, errors, t(Qs))
    )
    colnames(out) <- c(
      ".estimate",
      ".error",
      paste0(".q", 100 * probs)
    )
    class(out) <- c("tbl_df", "tbl", "data.frame")
  } else {
    out <- distributional::dist_sample(
      lapply(seq_len(NCOL(preds)), function(i) preds[, i])
    )
  }

  return(out)
}

#' Normalise quantiles into sampling weights
#' @noRd
norm_quantiles <- function(x) {
  xhat <- vector(length = length(x))
  for (i in seq_along(x)) {
    if (x[i] < 50) {
      xhat[i] <- (50 - x[i]) / 50
    } else {
      xhat[i] <- (x[i] - 50) / 50
    }
  }
  1.1 - xhat
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
