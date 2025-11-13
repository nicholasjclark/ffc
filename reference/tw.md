# Tweedie family

Tweedie distribution family for continuous data with zeros. Use in
[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)
`family` argument.

## Usage

``` r
tw(theta = NULL, link = "log", a = 1.01, b = 1.99)
```

## See also

[`tw`](https://rdrr.io/pkg/mgcv/man/Tweedie.html)

## Examples

``` r
# ffc_gam(continuous_with_zeros ~ fts(x), data = data, family = tw(), time = "time")
```
