# Beta regression family

Beta distribution family for data on (0,1) interval. Use in
[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)
`family` argument.

## Usage

``` r
betar(theta = NULL, link = "logit", eps = .Machine$double.eps * 100)
```

## See also

[`betar`](https://rdrr.io/pkg/mgcv/man/Beta.html)

## Examples

``` r
# ffc_gam(proportions ~ fts(x), data = data, family = betar(), time = "time")
```
