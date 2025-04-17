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
n_draws = 10
n_sims = 10
model = ETS()

# Extract the full linear predictor matrix
orig_lpmat <- predict(
  object,
  newdata = newdata,
  type = 'lpmatrix'
)

orig_betas <- mgcv::rmvn(
  n = n_draws * n_sims,
  mu = coef(object),
  V = vcov(object)
)

if(is.null(object$gam_init)) {

  # No need to modify lpmatrix if there were no
  # fts() terms in the model
  full_lpmat <- orig_lpmat
  full_betas <- orig_betas

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

  # Fit time series model to the basis coefficients
  # and generate forecasts
  functional_fc <- forecast(
    object = functional_coefs,
    h = max_horizon,
    times = n_sims,
    model = mod_name
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

  # Drop these columns from the lpmatrix
  intermed_lpmat <- orig_lpmat[, -fts_coefs, drop = FALSE]

  # Add columns of 1s for the fts basis coefficients
  full_lpmat <- cbind(
    intermed_lpmat,
    matrix(1,
           nrow = NROW(intermed_lpmat),
           ncol = NCOL(interpreted$data[grep('fts_',
                                             colnames(interpreted$data))])
    )
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
    )

  # Should now have n_draws * n_sims draws for each row of newdata
  if(!NROW(full_lpmat) * n_draws * n_sims == NROW(fts_fc)){
    stop('Wrong dimensions in forecast coefs; need to check on this')
  }
}
