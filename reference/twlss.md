# Tweedie location-scale-shape family

Extended Tweedie family allowing shape and scale parameters to vary. Use
in
[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)
`family` argument for complex Tweedie models.

## Usage

``` r
twlss(link = list("log", "identity", "identity"), a = 1.01, b = 1.99)
```

## See also

[`twlss`](https://rdrr.io/pkg/mgcv/man/twlss.html),
[`tw()`](https://nicholasjclark.github.io/ffc/reference/tw.md)

## Examples

``` r
# ffc_gam(list(y ~ fts(x), ~ s(z), ~ s(w)), family = twlss(), data = data, time = "time")
```
