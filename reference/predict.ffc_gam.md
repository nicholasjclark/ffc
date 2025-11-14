# Predict from a fitted ffc_gam model

Predict from a fitted ffc_gam model

## Usage

``` r
# S3 method for class 'ffc_gam'
predict(object, newdata, type = "link", se.fit = FALSE, ...)
```

## Arguments

- object:

  a fitted `gam` object as produced by
  [`gam()`](https://rdrr.io/pkg/mgcv/man/gam.html).

- newdata:

  A data frame or list containing the values of the model covariates at
  which predictions are required. If this is not provided then
  predictions corresponding to the original data are returned. If
  `newdata` is provided then it should contain all the variables needed
  for prediction: a warning is generated if not. See details for use
  with `link{linear.functional.terms}`.

- type:

  When this has the value `"link"` (default) the linear predictor
  (possibly with associated standard errors) is returned. When
  `type="terms"` each component of the linear predictor is returned
  seperately (possibly with standard errors): this includes parametric
  model components, followed by each smooth component, but excludes any
  offset and any intercept. `type="iterms"` is the same, except that any
  standard errors returned for smooth components will include the
  uncertainty about the intercept/overall mean. When `type="response"`
  predictions on the scale of the response are returned (possibly with
  approximate standard errors). When `type="lpmatrix"` then a matrix is
  returned which yields the values of the linear predictor (minus any
  offset) when postmultiplied by the parameter vector (in this case
  `se.fit` is ignored). The latter option is most useful for getting
  variance estimates for quantities derived from the model: for example
  integrated quantities, or derivatives of smooths. A linear predictor
  matrix can also be used to implement approximate prediction outside
  `R` (see example code, below).

- se.fit:

  when this is TRUE (not default) standard error estimates are returned
  for each prediction.

- ...:

  ignored

## Value

If `type == "lpmatrix"` then a `matrix` is returned which will give a
vector of linear predictor values (minus any offest) at the supplied
covariate values, when applied to the model coefficient vector.
Otherwise, if `se.fit == TRUE` then a `2 `item `list` is returned with
items (both `arrays`) `fit` and `se.fit` containing predictions and
associated standard error estimates, otherwise an `array` of predictions
is returned. The dimensions of the returned `arrays` depends on whether
type is `"terms"` or not: if it is then the `array` is `2` dimensional
with each term in the linear predictor separate, otherwise the `array`
is 1 dimensional and contains the linear predictor/predicted values (or
corresponding s.e.s). The linear predictor returned termwise will not
include the offset or the intercept.

## Details

This function returns predictions from models fitted with
[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md).
Data passed to `newdata` will first be correctly augmented to include
any basis functions whose coefficients were estimated as time-varying,
so the user need only supply data that includes variables that were used
in the original `data` that was supplied to
[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)

## See also

[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md),
[`fts()`](https://nicholasjclark.github.io/ffc/reference/fts.md)

## Author

Nicholas J Clark

## Examples

``` r
# Fit a model and generate predictions
mod <- ffc_gam(
  deaths ~ offset(log(population)) + sex + 
    fts(age, k = 8, bs = "cr", time_k = 10),
  time = "year", 
  data = qld_mortality,
  family = poisson(),
  engine = "bam"
)
preds <- predict(mod, type = "response")
```
