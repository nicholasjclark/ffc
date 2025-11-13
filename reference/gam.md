# Generalized Additive Model fitting

Fits GAMs using mgcv. In ffc, use
[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)
instead which extends gam functionality with time-varying coefficients.

## Usage

``` r
gam(
  formula,
  family = gaussian(),
  data = list(),
  weights = NULL,
  subset = NULL,
  na.action,
  offset = NULL,
  method = "GCV.Cp",
  optimizer = c("outer", "newton"),
  control = list(),
  scale = 0,
  select = FALSE,
  knots = NULL,
  sp = NULL,
  min.sp = NULL,
  H = NULL,
  gamma = 1,
  fit = TRUE,
  paraPen = NULL,
  G = NULL,
  in.out = NULL,
  drop.unused.levels = TRUE,
  drop.intercept = NULL,
  nei = NULL,
  discrete = FALSE,
  ...
)
```

## See also

[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html),
[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)

## Examples

``` r
# Use ffc_gam() for time-varying functionality
# mod <- ffc_gam(y ~ fts(x), data = data, time = "time")
```
