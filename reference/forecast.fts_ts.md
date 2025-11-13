# Forecasting functional basis coefficients

Forecasting functional basis coefficients

## Usage

``` r
# S3 method for class 'fts_ts'
forecast(object, model = "ARIMA", h = 1, times = 25, stationary = FALSE, ...)
```

## Arguments

- object:

  An object of class `fts_ts` containing time-varying basis function
  coefficients extracted from an `ffc_gam` object

- model:

  A `character` string representing a valid univariate model definition
  from the fable package or one of the built-in Bayesian dynamic factor
  models. Note that if a fable model is used, the chosen method must
  have an associated
  [`generate()`](https://nicholasjclark.github.io/ffc/reference/generate.md)
  method in order to simulate forecast realisations. Valid models
  currently include: `'ARDF'`, `'GPDF'`, '`VARDF`, `'ETS'`, `'ARIMA'`,
  `'AR'`, `'RW'`, `'NAIVE'`, and `'NNETAR'`

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
# \donttest{
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
forecast(coefs, h = 3, model = "ARIMA")
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> # A tsibble: 2,625 x 6 [1Y]
#> # Key:       .basis, .realisation, .model, .rep [875]
#>    .basis          .realisation .model .time .rep   .sim
#>    <chr>                  <int> <chr>  <dbl> <chr> <dbl>
#>  1 fts_bs_s_age__1            1 ARIMA   2021 1     -1.93
#>  2 fts_bs_s_age__1            1 ARIMA   2022 1     -1.87
#>  3 fts_bs_s_age__1            1 ARIMA   2023 1     -1.81
#>  4 fts_bs_s_age__1            1 ARIMA   2021 10    -1.93
#>  5 fts_bs_s_age__1            1 ARIMA   2022 10    -1.88
#>  6 fts_bs_s_age__1            1 ARIMA   2023 10    -1.83
#>  7 fts_bs_s_age__1            1 ARIMA   2021 11    -1.93
#>  8 fts_bs_s_age__1            1 ARIMA   2022 11    -1.87
#>  9 fts_bs_s_age__1            1 ARIMA   2023 11    -1.81
#> 10 fts_bs_s_age__1            1 ARIMA   2021 12    -1.93
#> # ℹ 2,615 more rows
# }
```
