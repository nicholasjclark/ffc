# Forecasting functional basis coefficients (Internal)

**For internal use:** This function is primarily used internally by
[`forecast.ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/forecast.ffc_gam.md).
Most users should call
[`forecast()`](https://generics.r-lib.org/reference/forecast.html)
directly on their `ffc_gam` object instead of using this function, which
requires properly structured coefficient data and has specific format
requirements.

## Usage

``` r
# S3 method for class 'fts_ts'
forecast(
  object,
  model = "ARIMA",
  h = get_stan_param("h", "forecast"),
  times = 25,
  stationary = FALSE,
  ...
)
```

## Arguments

- object:

  An object of class `fts_ts` containing time-varying basis function
  coefficients extracted from an `ffc_gam` object using
  [`fts_coefs()`](https://nicholasjclark.github.io/ffc/reference/fts_coefs.ffc_gam.md)

- model:

  A `character` string representing a valid univariate model definition
  from the fable package, ensemble methods, or one of the built-in
  Bayesian dynamic factor models. Note that if a fable model is used,
  the chosen method must have an associated
  [`generate()`](https://generics.r-lib.org/reference/generate.html)
  method in order to simulate forecast realisations. Valid models
  currently include: `'ENS'`, `'ARDF'`, `'GPDF'`, '`VARDF`, `'ETS'`,
  `'ARIMA'`, `'AR'`, `'RW'`, `'NAIVE'`, and `'NNETAR'`. The `'ENS'`
  option combines ETS and Random Walk forecasts with equal weights,
  hedging bets between exponential smoothing and random walk assumptions
  to provide more robust predictions when model uncertainty is high.

- h:

  A positive `integer` specifying the length of the forecast horizon

- times:

  A positive `integer` specifying the number of forecast realisation
  paths to simulate from the fitted forecast `model`

- stationary:

  If `TRUE`, the fitted time series models are constrained to be
  stationary. Default is `FALSE`. This option only works when
  `model == 'ARIMA'`

- ...:

  Other arguments to pass to the Stan dynamic factor models (i.e. the
  (V)AR order `lag = ...` or the number of factors `K = ...`)

## Value

A `tsibble` object containing the forecast prediction (`.sim`) for each
replicate realisation (`.rep`) at each timestep in the forecast horizon
`h`

## Details

Basis function coefficient time series will be used as input to the
specified `model` to train a forecasting model that will then be
extrapolated `h` timesteps into the future. A total of `times` forecast
realisations will be returned for each basis coefficient

## See also

[`predict()`](https://rdrr.io/r/stats/predict.html),
[`fts()`](https://nicholasjclark.github.io/ffc/reference/fts.md)

## Author

Nicholas J Clark

## Examples

``` r
# Extract coefficients and generate forecasts
mod <- ffc_gam(
  deaths ~ offset(log(population)) + sex +
    fts(age, k = 8, bs = "cr", time_k = 10),
  time = "year",
  data = qld_mortality,
  family = poisson(),
  engine = "bam"
)
coefs <- fts_coefs(mod, summary = FALSE, times = 5)

# Generate ETS forecasts
forecast(coefs, model = "ETS", h = 3)
#> # A tsibble: 2,625 x 6 [1Y]
#> # Key:       .basis, .realisation, .model, .rep [875]
#>    .basis          .realisation .model  year .rep   .sim
#>    <chr>                  <int> <chr>  <dbl> <chr> <dbl>
#>  1 fts_bs_s_age__1            1 ETS     2021 1     -1.90
#>  2 fts_bs_s_age__1            1 ETS     2022 1     -1.85
#>  3 fts_bs_s_age__1            1 ETS     2023 1     -1.79
#>  4 fts_bs_s_age__1            1 ETS     2021 10    -1.91
#>  5 fts_bs_s_age__1            1 ETS     2022 10    -1.86
#>  6 fts_bs_s_age__1            1 ETS     2023 10    -1.81
#>  7 fts_bs_s_age__1            1 ETS     2021 11    -1.90
#>  8 fts_bs_s_age__1            1 ETS     2022 11    -1.84
#>  9 fts_bs_s_age__1            1 ETS     2023 11    -1.77
#> 10 fts_bs_s_age__1            1 ETS     2021 12    -1.91
#> # ℹ 2,615 more rows
```
