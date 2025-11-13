# Negative binomial family (alternative)

Alternative negative binomial family specification. Use in
[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)
`family` argument.

## Usage

``` r
negbin(theta = stop("'theta' must be specified"), link = "log")
```

## See also

[`negbin`](https://rdrr.io/pkg/mgcv/man/negbin.html)

## Examples

``` r
# ffc_gam(counts ~ fts(x), data = data, family = negbin(), time = "time")
```
