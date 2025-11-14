# Forecasting `ffc_gam` models

Forecasting `ffc_gam` models

## Usage

``` r
# S3 method for class 'ffc_gam'
forecast(
  object,
  newdata,
  type = "response",
  model = "ARIMA",
  stationary = FALSE,
  summary = TRUE,
  robust = TRUE,
  probs = c(0.025, 0.1, 0.9, 0.975),
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

- ...:

  Other arguments to pass to the Stan dynamic factor models

## Value

Predicted values on the appropriate scale. If `summary == FALSE`, the
output is a matrix. If `summary == TRUE`, the output is a tidy
`tbl_df / data.frame`

## Details

Computes forecast distributions from fitted `ffc_gam` objects

## See also

[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md),
[`fts()`](https://nicholasjclark.github.io/ffc/reference/fts.md),
[`forecast.fts_ts()`](https://nicholasjclark.github.io/ffc/reference/forecast.fts_ts.md)

## Author

Nicholas J Clark
