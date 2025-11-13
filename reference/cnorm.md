# Censored normal family

Censored normal distribution family for censored continuous data. Use in
[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)
`family` argument.

## Usage

``` r
cnorm(theta = NULL, link = "identity")
```

## See also

[`cnorm`](https://rdrr.io/pkg/mgcv/man/cnorm.html)

## Examples

``` r
# ffc_gam(censored_y ~ fts(x), data = data, family = cnorm(), time = "time")
```
