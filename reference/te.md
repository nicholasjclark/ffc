# Tensor product smooth terms

Constructs tensor product smooth terms for multidimensional smoothing.
Use within
[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)
formulae for interactions between variables.

## Usage

``` r
te(
  ...,
  k = NA,
  bs = "cr",
  m = NA,
  d = NA,
  by = NA,
  fx = FALSE,
  np = TRUE,
  xt = NULL,
  id = NULL,
  sp = NULL,
  pc = NULL
)
```

## See also

[`te`](https://rdrr.io/pkg/mgcv/man/te.html)

## Examples

``` r
# Tensor product smooth in ffc_gam
# ffc_gam(y ~ te(x1, x2), data = data, time = "time")
```
