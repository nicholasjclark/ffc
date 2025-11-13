# Tweedie log density

Computes log density for Tweedie distribution. Utility function for
Tweedie family.

## Usage

``` r
ldTweedie(
  y,
  mu = y,
  p = 1.5,
  phi = 1,
  rho = NA,
  theta = NA,
  a = 1.001,
  b = 1.999,
  all.derivs = FALSE
)
```

## See also

[`ldTweedie`](https://rdrr.io/pkg/mgcv/man/ldTweedie.html),
[`tw()`](https://nicholasjclark.github.io/ffc/reference/tw.md)

## Examples

``` r
# ldTweedie(y = data, mu = fitted_values, p = 1.5, phi = 1)
```
