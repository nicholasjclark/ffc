# Forecasting `ffc_gam` models

Forecasting `ffc_gam` models

## Usage

``` r
# S3 method for class 'ffc_gam'
forecast(
  object,
  newdata,
  type = "response",
  model = "ENS",
  stationary = FALSE,
  summary = TRUE,
  robust = TRUE,
  probs = c(0.025, 0.1, 0.9, 0.975),
  mean_model = "ENS",
  ...
)
```

## Arguments

- object:

  An object of class `ffc_gam`. See
  [`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)

- newdata:

  `dataframe` or `list` of test data containing the same variables that
  were included in the original `data` used to fit the model. The
  covariate information in `newdata`, along with the temporal
  information indexed by the `time` variable in the original call to
  [`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md),
  will be used to generate forecasts from the fitted model equations

- type:

  When this has the value `link`, the linear predictor is calculated on
  the link scale. If `expected` is used, predictions reflect the
  expectation of the response (the mean) but ignore uncertainty in the
  observation process. When `response` is used (the default), the
  predictions take uncertainty in the observation process into account
  to return predictions on the outcome scale

- model:

  A character string specifying the forecasting model to use. Default is
  "ENS" (ensemble). Options include "ENS", "ETS", "ARIMA", "RW",
  "NAIVE", and Stan dynamic factor models ("ARDF", "VARDF", "GPDF").
  "ENS" combines ETS and Random Walk forecasts with equal weights,
  hedging bets between different forecasting assumptions to improve
  robustness.

- stationary:

  If `TRUE`, the fitted time series models are constrained to be
  stationary. Default is `FALSE`. This option only works when
  `model == 'ARIMA'`

- summary:

  Should summary statistics be returned instead of the raw values?
  Default is `TRUE`

- robust:

  If `FALSE` (the default) the mean is used as the measure of central
  tendency and the standard deviation as the measure of variability. If
  `TRUE`, the median and the median absolute deviation (MAD) are applied
  instead. Only used if `summary` is `TRUE`

- probs:

  The percentiles to be computed by the
  [`quantile()`](https://rdrr.io/r/stats/quantile.html) function. Only
  used if `summary` is `TRUE`

- mean_model:

  A character string specifying the forecasting model to use for mean
  basis coefficients when using Stan factor models (ARDF, VARDF, GPDF).
  Default is "ENS". Options include "ENS", "ETS", "ARIMA", "RW",
  "NAIVE". "ENS" creates an ensemble of ETS and RW forecasts with equal
  weights, hedging bets between different forecasting assumptions for
  mean coefficients. This is only used when forecasting mixed
  mean/non-mean basis functions with Stan factor models.

- ...:

  Additional arguments for Stan dynamic factor models (ARDF, VARDF,
  GPDF). Key arguments include: K (number of factors, default: 2), lag
  (AR order, default: 1), chains (MCMC chains, default: 4), iter
  (iterations, default: 500), silent (suppress progress, default: TRUE),
  cores, adapt_delta, max_treedepth.

## Value

Predicted values on the appropriate scale. If `summary == FALSE`, the
output is a matrix. If `summary == TRUE`, the output is a tidy
`tbl_df / data.frame`

## Details

**Distributional Family Support:**

Full support for distributional regression models using mgcv families:

- **gaulss:** Gaussian location-scale (normal distribution with varying
  mean and variance)

- **twlss:** Tweedie location-scale-shape

For distributional families, forecasting operates on all parameters
simultaneously, producing forecasts that capture both mean and variance
dynamics over time.

**Forecasting Methodology:**

This function implements a two-stage forecasting approach for functional
regression models with time-varying coefficients of the form: \$\$y_t =
\sum\_{j=1}^{J} \beta_j(t) B_j(x) + \epsilon_t\$\$ where \\\beta_j(t)\\
are time-varying coefficients and \\B_j(x)\\ are basis functions.

1.  **Extract basis coefficients:** Time-varying functional coefficients
    \\\beta_j(t)\\ are extracted from the fitted GAM as time series

2.  **Forecast coefficients:** These coefficient time series are
    forecast using ensemble methods (ENS), Stan dynamic factor models
    (ARDF/VARDF/GPDF), or individual fable models (ETS, ARIMA, etc.)

3.  **Reconstruct forecasts:** Forecasted coefficients are combined:
    \$\$\hat{y}\_{t+h} = \sum\_{j=1}^{J} \hat{\beta}\_j(t+h) B_j(x)\$\$

4.  **Combine uncertainties:** Multiple uncertainty sources are
    integrated hierarchically (see Uncertainty Quantification section)

**Uncertainty Quantification:**

Forecast uncertainty is captured through a hierarchical structure:

*Within Stan dynamic factor models:*

- **Process uncertainty:** Factor dynamics, autoregressive terms, factor
  loadings

- **Observation uncertainty:** Series-specific error terms
  (\\\sigma\_{obs}\\)

*Final combination in linear predictor space:*

- **Stan forecast samples:** Already incorporate process + observation
  uncertainty

- **GAM parameter uncertainty:** Random draws from
  \\N(\hat{\boldsymbol{\theta}}, \mathbf{V})\\ where \\\mathbf{V}\\ is
  the coefficient covariance matrix

These components are combined additively: \\\text{Stan forecasts} +
\text{GAM uncertainty}\\

**Model Selection:**

- **ENS ensemble (default):** Combines ETS and Random Walk forecasts
  with equal weights. This hedges bets between exponential smoothing
  assumptions (trend and seasonality patterns continue) and random walk
  assumptions (future values equal current values). Provides robust
  predictions when model uncertainty is high, which is common in
  coefficient forecasting.

- **Stan factor models (ARDF/VARDF/GPDF):** Used for multivariate
  forecasting of non-mean basis coefficients. Capture dependencies
  between coefficient series and assume zero-centered time series for
  efficiency.

- **Mean basis models:** Used for mean basis coefficients (which operate
  at non-zero levels). Default is ENS ensemble, controlled by
  `mean_model` parameter.

**Important Note on `times` Parameter:**

For Stan dynamic factor models, the `times` parameter is automatically
set to `(iter - warmup) * chains` to ensure dimensional consistency. Any
user-specified `times` value will be ignored with a warning. For ARIMA
models, `times` can be specified freely and controls the number of
posterior draws.

## See also

[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md),
[`fts()`](https://nicholasjclark.github.io/ffc/reference/fts.md),
[`forecast.fts_ts()`](https://nicholasjclark.github.io/ffc/reference/forecast.fts_ts.md)

## Author

Nicholas J Clark

## Examples

``` r
# Basic forecasting example with growth data
data("growth_data")
mod <- ffc_gam(
  height_cm ~ s(id, bs = "re") +
    fts(age_yr, k = 8, bs = "cr", time_k = 10),
  time = "age_yr", 
  data = growth_data,
  family = Gamma()
)

# Forecast with ETS model
newdata <- data.frame(
  id = "boy_11", 
  age_yr = c(16, 17, 18)
)
fc <- forecast(mod, newdata = newdata)  # Uses ENS ensemble by default

# Forecast with specific models
fc_ets <- forecast(mod, newdata = newdata, model = "ETS")
fc_rw <- forecast(mod, newdata = newdata, model = "RW")

# Get raw forecast matrix without summary
fc_raw <- forecast(mod, newdata = newdata, summary = FALSE)
```
