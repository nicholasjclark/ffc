#' Fit an autoregressive dymamic factor model
#'
#' Fit an autoregressive dymamic factor model using Stan
#'
#' @param data a `list` of data to pass to Stan
#' @param chains `integer` specifying the number of chains to be run
#' @param iter `integer` specifying the number of iterations per chain
#' @param warmup `integer` specifying the number of initial iterations to
#' use as burnin
#' @param adapt_delta the thin of the jumps in a HMC method
#' @param max_treedepth maximum tree depth per iteration
#'
#' @author Nicholas J Clark
#'
#' @noRd
#'
#' @return a stanfit object
#'
fit_ardf = function(data,
                    chains = 4,
                    iter = 1000,
                    warmup = floor(iter / 2),
                    adapt_delta = 0.80,
                    max_treedepth = 9,
                    ...){

  stanfit <- rstan::sampling(
    stanmodels$ardf,
    data = data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    control = list(adapt_delta = adapt_delta,
                   max_treedepth = max_treedepth)
  )

  return(stanfit)
}
