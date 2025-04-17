# Forecasting from a ffc_gam object
qld_train <- qld_mortality %>%
  dplyr::filter(year < 2016)
qld_test <- qld_mortality %>%
  dplyr::filter(year >= 2016)
mod <- ffc_gam(
  deaths ~
    offset(log(population)) +
    sex +
    fts(
      age,
      k = 10, bs = "cr", by = sex,
      time_bs = "cr", time_k = 10
    ),
  time = "year",
  data = qld_train,
  family = poisson(),
  engine = "bam"
)
summary(mod)

functional_fc <- forecast(
  object = functional_coefs,
  h = 5,
  times = 5,
  model = 'ETS'
)
functional_fc


object = mod
newdata = qld_test
n_draws = 5
n_sims = 10
model = ETS()

# Extract the full linear predictor matrix
orig_lpmat <- predict(
  object,
  newdata = newdata,
  type = 'lpmatrix'
)

# Take full draws of beta coefficients
orig_betas <- mgcv::rmvn(
  n = n_draws * n_sims,
  mu = coef(object),
  V = vcov(object)
)

if(is.null(object$gam_init)) {

  # No need to modify lpmatrix if there were no
  # fts() terms in the model
  full_linpreds <- matrix(
    as.vector(
      t(apply(
        as.matrix(orig_betas),
        1,
        function(row) orig_lpmat %*% row +
          attr(orig_lpmat, 'model.offset')
      ))
    ),
    nrow = NROW(orig_betas)
  )

} else {
  # Determine horizon (assuming equal time gaps)
  time_var <- object$time_var

  interpreted <- ffc:::interpret_ffc(
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
  functional_coefs <- fts_coefs(
    object,
    summary = FALSE,
    times = n_draws
  )

  # Validate model choice
  if(is.character(model)) {
    mod_name <- model
  } else {
    mod_name <- deparse(substitute(model))
    mod_name <- gsub(
      "\\s*(\\([^()]*(?:(?1)[^()]*)*\\))",
      "",
      mod_name,
      perl = TRUE
    )
  }

  # Fit the time series model to the basis coefficients
  # and generate forecasts
  functional_fc <- forecast(
    object = functional_coefs,
    h = max_horizon,
    times = n_sims,
    model = mod_name
  )

  # Only need to return forecasts for those times that are
  # in newdata
  functional_fc <- functional_fc %>%
    dplyr::filter(
      .time %in% unique(interpreted$data[[time_var]])
    )

  # Which coefficients in lpmatrix are associated with fts objects?
  smooth_names <- unlist(
    purrr::map(object$smooth, 'label'),
    use.names = FALSE
  )

  fts_names <- grep(
    ':fts_',
    smooth_names,
    fixed = TRUE
  )

  fts_coefs <- unlist(
    purrr::map(
      object$smooth[fts_names],
      \ (x) x$first.para : x$last.para
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
        function(row) intermed_lpmat %*% row +
          attr(orig_lpmat, 'model.offset')
      ))
    ),
    nrow = NROW(intermed_betas)
  )

  # Join forecasts to the basis function evaluations
  fts_fc <- data.frame(
    interpreted$data[grep('fts_', colnames(interpreted$data))]
  ) %>%
    dplyr::bind_cols(
      data.frame(.time = interpreted$data[[time_var]],
                 .row = 1:length(interpreted$data[[1]]))
    ) %>%
    tidyr::pivot_longer(
      cols = !contains(c(".time",
                         ".row")),
      names_to = ".basis",
      values_to = '.evaluation'
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
    dplyr::mutate(.draw = paste0(.realisation, '_', .rep)) %>%
    dplyr::select(-.time,
                  -.realisation,
                  -.rep)

  # Should now have n_draws * n_sims draws for each row of newdata
  if(!NROW(full_lpmat) * n_draws * n_sims == NROW(fts_fc)){
    stop('Wrong dimensions in forecast coefs; need to check on this')
  }

  # Should also have same ncols as number of fts basis functions
  if(
    !NCOL(fts_fc) - 2 ==
    NCOL(interpreted$data[
      grep('fts_',
           colnames(interpreted$data))])
  ){
    stop('Wrong dimensions in forecast coefs; need to check on this')
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
      function(x){
        fc_linpreds[which(fts_fc$.draw == x)]
      }
    )
  ) +
    intermed_linpreds
}

# Now can proceed to send full_linpreds to the relevant
# invlink and rng functions for outcome-level predictions;
# use gratia::posterior_samples() codebase as a guide
