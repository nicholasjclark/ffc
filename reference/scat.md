# Scaled Student-t family

Scaled Student-t distribution family for robust regression. Use in
[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)
`family` argument or Stan models.

## Usage

``` r
scat(theta = NULL, link = "identity", min.df = 3)
```

## See also

[`scat`](https://rdrr.io/pkg/mgcv/man/scat.html)

## Examples

``` r
# ffc_gam(y ~ fts(x), data = data, family = scat(), time = "time")
```
