# Plot time series of fts basis coefficients

Produces a time series plot of basis function coefficients from `fts_ts`
objects

## Usage

``` r
# S3 method for class 'fts_ts'
autoplot(object, ...)
```

## Arguments

- object:

  An object of class `fts_ts` containing time-varying basis function
  coefficients extracted from an `ffc_gam` object

- ...:

  Ignored

## Author

Nicholas J Clark

## Examples

``` r
# Extract and plot time-varying coefficients
mod <- ffc_gam(
  deaths ~ offset(log(population)) + sex + 
    fts(age, k = 8, bs = "cr", time_k = 10),
  time = "year", 
  data = qld_mortality,
  family = poisson(),
  engine = "bam"
)
coefs <- fts_coefs(mod, summary = FALSE, n_samples = 5)
autoplot(coefs)
```
