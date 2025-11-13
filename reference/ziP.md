# Zero-inflated Poisson family

Zero-inflated Poisson distribution family for count data with excess
zeros. Use in
[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)
`family` argument.

## Usage

``` r
ziP(theta = NULL, link = "identity", b = 0)
```

## See also

[`ziP`](https://rdrr.io/pkg/mgcv/man/ziP.html)

## Examples

``` r
# ffc_gam(zero_inflated_counts ~ fts(x), data = data, family = ziP(), time = "time")
```
