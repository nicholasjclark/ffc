# Extract time-varying basis coefficients

Extract time-varying basis coefficients a fitted ffc_gam model so they
can be visualized and / or forecasted ahead

## Usage

``` r
fts_coefs(object, ...)

# S3 method for class 'ffc_gam'
fts_coefs(object, summary = TRUE, times = 25, ...)
```

## Arguments

- object:

  `list` object of class `ffc_gam`. See
  [`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)

- ...:

  Ignored

- summary:

  `Logical`. Should summary statistics of the coefficient time series be
  returned instead of realized curves? Default is `TRUE`. If `FALSE`,
  replicate realisations of each basis coefficient time series will be
  returned

- times:

  A positive `integer` specifying the number of time series realisation
  paths to simulate from the fitted model. Ignored if `summary = FALSE`

## Value

A `fts_ts` object containing the point estimates and their standard
errors for basis function coefficients

## Details

This function creates a `tidy` time series of basis function
coefficients from all
[`fts()`](https://nicholasjclark.github.io/ffc/reference/fts.md) terms
that were supplied to the original model

## See also

[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md),
[`fts()`](https://nicholasjclark.github.io/ffc/reference/fts.md)

## Author

Nicholas J Clark
