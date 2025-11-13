# Ordered categorical family

Ordered categorical distribution family for ordinal response data. Use
in
[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)
`family` argument.

## Usage

``` r
ocat(theta = NULL, link = "identity", R = NULL)
```

## See also

[`ocat`](https://rdrr.io/pkg/mgcv/man/ocat.html)

## Examples

``` r
# ffc_gam(ordered_response ~ fts(x), data = data, family = ocat(), time = "time")
```
