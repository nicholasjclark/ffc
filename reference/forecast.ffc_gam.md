# Forecasting `ffc_gam` models

Forecasting `ffc_gam` models

## Usage

``` r
# S3 method for class 'ffc_gam'
forecast(
  object,
  newdata,
  type = "response",
  model = "ETS",
  stationary = FALSE,
  summary = TRUE,
  robust = TRUE,
  probs = c(0.025, 0.1, 0.9, 0.975),
  mean_model = "ETS",
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

  A `character` string representing a valid univariate model definition
  from the fable package or one of the built-in Bayesian dynamic factor
  models. Note that if a fable model is used, the chosen method must
  have an associated
  [`generate()`](https://generics.r-lib.org/reference/generate.html)
  method in order to simulate forecast realisations. Valid models
  currently include: `'ARDF'`, `'GPDF'`, '`VARDF`, `'ETS'`, `'ARIMA'`,
  `'AR'`, `'RW'`, `'NAIVE'`, and `'NNETAR'`

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
  Default is "ETS". Options include "ETS", "ARIMA", "RW", "NAIVE". This
  is only used when forecasting mixed mean/non-mean basis functions with
  Stan factor models.

- ...:

  Additional arguments for Stan dynamic factor models (ARDF, VARDF,
  GPDF):

  K

  :   Number of latent factors (default: 2). Must be positive.

  lag

  :   Number of time lags for autoregressive terms (default: 1). Must be
      positive.

  chains

  :   Number of MCMC chains (default: 4)

  iter

  :   Total iterations per chain (default: 500)

  warmup

  :   Warmup iterations per chain (default: iter/2)

  cores

  :   Number of CPU cores to use (default: min(chains, available cores))

  adapt_delta

  :   Target acceptance rate (default: 0.75)

  max_treedepth

  :   Maximum tree depth (default: 9)

  times

  :   Number of posterior samples to draw for coefficient forecasting.
      For Stan models, this is automatically set to (iter - warmup) \*
      chains to ensure dimension consistency. Minimum 100 total
      posterior samples required.

## Value

Predicted values on the appropriate scale. If `summary == FALSE`, the
output is a matrix. If `summary == TRUE`, the output is a tidy
`tbl_df / data.frame`

## Details

**Forecasting Methodology:**

This function implements a two-stage forecasting approach for functional
regression models with time-varying coefficients of the form: \$\$y_t =
\sum\_{j=1}^{J} \beta_j(t) B_j(x) + \epsilon_t\$\$ where \\\beta_j(t)\\
are time-varying coefficients and \\B_j(x)\\ are basis functions.

1.  **Extract basis coefficients:** Time-varying functional coefficients
    \\\beta_j(t)\\ are extracted from the fitted GAM as time series

2.  **Forecast coefficients:** These coefficient time series are
    forecast using either Stan dynamic factor models (ARDF/VARDF/GPDF)
    or ARIMA models

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

- **Stan factor models (ARDF/VARDF/GPDF):** Used for multivariate
  forecasting of non-mean basis coefficients. Capture dependencies
  between coefficient series and assume zero-centered time series for
  efficiency.

- **Mean basis models:** Used for mean basis coefficients (which operate
  at non-zero levels). Default is ETS, controlled by `mean_model`
  parameter.

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
# \donttest{
# Basic forecasting example
data("growth_data")
mod <- ffc_gam(
  Reaction ~ fts(Subject, k = 3, time_k = 8) +
             fts(Days, k = 3, time_k = 8),
  time = "Days", data = growth_data
)
#> Error: the variable 'Days' cannot be found in data

# Forecast with ETS (default for mixed basis)
newdata <- data.frame(Subject = rep("308", 3), Days = 10:12)
fc <- forecast(mod, newdata = newdata, model = "ETS")
#> Error: object 'mod' not found

# Use ARDF with custom mean_model for mean basis
fc_ardf <- forecast(mod, newdata = newdata,
                    model = "ARDF",
                    mean_model = "RW",  # Use random walk for mean basis
                    K = 2)
#> Error: object 'mod' not found

# Get raw forecast matrix without summary
fc_raw <- forecast(mod, newdata = newdata,
                  model = "ETS",
                  summary = FALSE)
#> Error: object 'mod' not found
# }
```
