#' Create Stan data list from an ffc tbl_ts object
#' @noRd
prep_tbl_ts_stan = function(.data,
                            h,
                            K,
                            p,
                            family,
                            model){

  model <- match.arg(model, choices = c('ardf', 'gpdf', 'vardf'))

  # Expand the data to include NAs for the out of sample period
  .fulldata <- expand_tbl_ts(.data, h) %>%

    # add 'series' indicator
    dplyr::mutate(
      .series = factor(
        .basis,
        levels = set_series_levels(unique(.basis)
        )
      ),
      .y = .estimate
    ) %>%
    dplyr::arrange(.series, .time) %>%
    dplyr::select(.series, .time, .y, .sd)

  # Pull some key arguments
  family <- ifelse(family$family == 'gaussian', 1, 2)
  n_series <- nlevels(.fulldata$.series)

  if(K > n_series) {
    rlang::warn('K cannot be greater than the number of unique series. Setting K = n_series')
    K <- n_series
  }

  n_timepoints <- length(unique(.fulldata$.time))

  # Extract observation error for the last training time period
  # per series
  sample_sd <- .fulldata %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::group_by(.series) %>%
    dplyr::arrange(.time) %>%
    dplyr::summarise(.finalsd = tail(.sd, 1)) %>%
    dplyr::pull(.finalsd)

  # Need a matrix of values for the series, excluding forecast values
  Y <- .fulldata %>%
    dplyr::select(-.sd) %>%
    tidyr::pivot_wider(names_from = '.series',
                       values_from = '.y') %>%
    dplyr::select(-.time) %>%
    as.matrix()

  stopifnot(identical(NROW(Y), n_timepoints))
  stopifnot(identical(NCOL(Y), n_series))

  # Calculate series means
  alpha <- apply(Y, 2, function(x) mean(x, na.rm = TRUE))

  # Create matrix representing whether an observation was missing or not
  Y_observed <- matrix(
    NA,
    ncol = NCOL(Y),
    nrow = NROW(Y)
  )
  for (i in 1:dim(Y)[1]) {
    for (s in 1:dim(Y)[2]) {
      if (is.na(Y[i, s])) {
        Y_observed[i, s] = 0
      } else {
        Y_observed[i, s] = 1
      }
    }
  }

  # Use -1 for any missing observations so Stan doesn't throw errors due to NAs
  Y[is.na(Y)] <- -1

  # Begin the Stan model data list
  model_data <- list(
    n = n_timepoints,
    K = K,
    n_series = n_series,
    sample_sd = sample_sd,
    M = K * (n_series - K) + K * (K - 1) / 2 + K,
    n_nonmissing = length(which(Y_observed == 1)),
    flat_ys = as.vector(Y)[which(
      as.vector(Y_observed) == 1
    )],
    obs_ind = which(as.vector(Y_observed) == 1),
    family = family,
    alpha = alpha,
    beta = array(1)
  )

  if (model %in% c('ardf', 'vardf')){
    if(p > max(.fulldata$.time) - 1){
      stop('Cannot estimate a model where p > maximum lag available in the data',
           call. = FALSE)
    }

    if (model == 'ardf'){
      c(model_data) <- nlist(
        prior_ar = c(1, 1, 1),
        P = p
      )
    }

    if (model == 'vardf'){
      c(model_data) <- nlist(
        P = p
      )
    }
  }

  return(model_data)
}

#' @noRd
extract_stan_fc = function(stanfit,
                           mod_name = 'ARDF',
                           .data,
                           model_data,
                           h){

  # Get forecast times
  min_time <- min(.data$.time)
  max_time <- max(.data$.time) + h

  # Extract forecast distributions for each series
  preds <- as.matrix(stanfit, pars = 'ypred')
  n_series <- model_data$n_series
  ends <- seq(
    0,
    NCOL(preds),
    length.out = n_series + 1
  )
  starts <- ends + 1
  starts <- c(1, starts[-c(1, (n_series + 1))])
  ends <- ends[-1]

  series_fcs <- do.call(
    rbind,
    lapply(
      1:n_series,
      function(x){
        s_name <- set_series_levels(unique(.data$.basis))[x]
        data.frame(
          .basis = s_name,
          .realisation = 1,
          .model = mod_name,
          .rep = 1 : NROW(preds),
          .sim = as.vector(preds[ , starts[x] : ends[x]][, -c(1 : (ends[1] - h))]),
          .time = sort(rep((ends[1] - h + 1) : (ends[1]),
                           NROW(preds)))
        )
      }
    )
  ) %>%
    dplyr::as_tibble()

  # Extend any other time variables ...
  index <- tsibble::index_var(.data)
  out <- suppressMessages(
    make_future_data(.data, h = h) %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(.time = rep((ends[1] - h + 1) : (ends[1]),
                                n_series)) %>%
      dplyr::select(-.basis, -.realisation) %>%
      dplyr::distinct() %>%

      # ... and join to the forecast data
      dplyr::right_join(series_fcs,
                        relationship = 'many-to-many',
                        by = dplyr::join_by(.time)) %>%
      tsibble::as_tsibble(
        key = c(.basis, .realisation, .model, .rep),
        index = index
      )
  )

  return(out)
}

#' @noRd
set_series_levels = function(x){
  unique_levels <- sort(unique(x))
  if(any(grepl('_mean', unique_levels))){
    newlevels <- vector(length = length(unique_levels))
    n_means <- length(which(grepl('_mean', unique_levels)))
    newlevels <- c(
      grep('_mean', unique_levels, value = TRUE),
      unique_levels[-grep('_mean', unique_levels)]
    )
  } else {
    newlevels <- unique_levels
  }
  return(newlevels)
}
