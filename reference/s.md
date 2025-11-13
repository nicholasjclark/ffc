# Smooth terms for use in GAMs

Constructs smooth terms for inclusion in GAM formulae. Use within
[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)
for standard smooths or
[`fts()`](https://nicholasjclark.github.io/ffc/reference/fts.md) for
time-varying smooths.

## Usage

``` r
s(
  ...,
  k = -1,
  fx = FALSE,
  bs = "tp",
  m = NA,
  by = NA,
  xt = NULL,
  id = NULL,
  sp = NULL,
  pc = NULL
)
```

## See also

[`s`](https://rdrr.io/pkg/mgcv/man/s.html)

## Examples

``` r
# Standard smooth in ffc_gam
# ffc_gam(y ~ s(x), data = data, time = "time")
```
