# Negative binomial family

Negative binomial distribution family for count data with
overdispersion. Use in
[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)
`family` argument.

## Usage

``` r
nb(theta = NULL, link = "log")
```

## See also

[`nb`](https://rdrr.io/pkg/mgcv/man/negbin.html)

## Examples

``` r
# ffc_gam(counts ~ fts(x), data = data, family = nb(), time = "time")
```
