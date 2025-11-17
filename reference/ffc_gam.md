# Fit a functional time series model using dynamic functional coefficients

Fit Generalized Additive Models that can include time-varying (dynamic)
functions

## Usage

``` r
ffc_gam(
  formula,
  family = gaussian(),
  data = list(),
  time,
  engine = c("gam", "bam"),
  ...
)
```

## Arguments

- formula:

  A GAM formula (see
  [`formula.gam`](https://rdrr.io/pkg/mgcv/man/formula.gam.html) and
  also [`gam.models`](https://rdrr.io/pkg/mgcv/man/gam.models.html)).
  This is exactly like the formula for a GLM except that smooth terms,
  [`fts()`](https://nicholasjclark.github.io/ffc/reference/fts.md),
  [`s()`](https://rdrr.io/pkg/mgcv/man/s.html) and
  [`te()`](https://rdrr.io/pkg/mgcv/man/te.html) can be added to the
  right hand side to specify that the linear predictor depends on smooth
  functions of predictors (or linear functionals of these).

- family:

  This is a family object specifying the distribution and link to use in
  fitting etc. See [`glm`](https://rdrr.io/r/stats/glm.html) and
  [`family`](https://rdrr.io/r/stats/family.html) for more details. The
  extended families listed in
  [`family.mgcv`](https://rdrr.io/pkg/mgcv/man/family.mgcv.html) can
  also be used.

- data:

  A data frame or list containing the model response variable and
  covariates required by the formula. By default the variables are taken
  from `environment(formula)`: typically the environment from which
  `gam` is called.

- time:

  `character` specifying which variable in `data` represents the the
  time ordering of the observations

- engine:

  `character` string specifying which mgcv interface to use for fitting
  the model.

- ...:

  other arguments to pass to either
  [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html)

## Value

An object of class `ffc_gam`, which inherits from objects of class `gam`
or `bam`

## Details

This function will update the supplied `formula` to ensure any
time-varying functionals (supplied through
[`fts()`](https://nicholasjclark.github.io/ffc/reference/fts.md) terms
in the formula right hand side) are appropriately incorporated into the
model. It then passes the updated model and data objects to the
specified `engine` for model fitting

## See also

[`fts()`](https://nicholasjclark.github.io/ffc/reference/fts.md),
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html),
[`bam`](https://rdrr.io/pkg/mgcv/man/bam.html)

## Author

Nicholas J Clark

## Examples

``` r
# Fit a dynamic function-on-function regression to the Queensland
# mortality data
data("qld_mortality")
mod <- ffc_gam(
  deaths ~
    offset(log(population)) +
    sex +
    fts(age,
      k = 8, bs = "cr",
      time_bs = "cr", time_k = 10
    ),
  time = "year",
  data = qld_mortality,
  family = poisson(),
  engine = "bam"
)
class(mod)
#> [1] "ffc_gam" "bam"     "gam"     "glm"     "lm"     
summary(mod)
#> 
#> Family: poisson 
#> Link function: log 
#> 
#> Formula:
#> deaths ~ sex + offset(log(population)) + s(year, by = fts_bs_s_age__1, 
#>     bs = "cr", k = 10, m = 2, id = 1) + s(year, by = fts_bs_s_age__2, 
#>     bs = "cr", k = 10, m = 2, id = 1) + s(year, by = fts_bs_s_age__3, 
#>     bs = "cr", k = 10, m = 2, id = 1) + s(year, by = fts_bs_s_age__4, 
#>     bs = "cr", k = 10, m = 2, id = 1) + s(year, by = fts_bs_s_age__5, 
#>     bs = "cr", k = 10, m = 2, id = 1) + s(year, by = fts_bs_s_age__6, 
#>     bs = "cr", k = 10, m = 2, id = 1) + s(year, by = fts_bs_s_age__7, 
#>     bs = "cr", k = 10, m = 2, id = 1)
#> 
#> Parametric coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -5.560332   0.002416 -2301.4   <2e-16 ***
#> sexmale      0.472635   0.002077   227.5   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Approximate significance of smooth terms:
#>                           edf Ref.df Chi.sq p-value    
#> s(year):fts_bs_s_age__1 5.546  6.517  62927  <2e-16 ***
#> s(year):fts_bs_s_age__2 6.157  7.219  47283  <2e-16 ***
#> s(year):fts_bs_s_age__3 6.547  7.634  26662  <2e-16 ***
#> s(year):fts_bs_s_age__4 7.460  8.526  10047  <2e-16 ***
#> s(year):fts_bs_s_age__5 8.158  9.096 196408  <2e-16 ***
#> s(year):fts_bs_s_age__6 8.152  9.046 642373  <2e-16 ***
#> s(year):fts_bs_s_age__7 6.170  7.233 639228  <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> R-sq.(adj) =  0.973   Deviance explained = 94.8%
#> fREML =  47425  Scale est. = 1         n = 8282

# Predictions work in the usual way
head(predict(mod, type = "link"))
#>        1        2        3        4        5        6 
#> 3.985387 4.010376 4.023095 3.998402 3.945344 3.902329 
head(predict(mod, type = "response"))
#>        1        2        3        4        5        6 
#> 53.80613 55.16762 55.87379 54.51097 51.69412 49.51764 

# Extract basis coefficient time series
(functional_ts <- fts_coefs(mod))
#> # A tibble: 287 × 5
#>    .basis          .time .estimate     .se  year
#>  * <chr>           <int>     <dbl>   <dbl> <int>
#>  1 fts_bs_s_age__1  1980     -4.13 0.00635  1980
#>  2 fts_bs_s_age__1  1981     -4.07 0.00521  1981
#>  3 fts_bs_s_age__1  1982     -4.02 0.00477  1982
#>  4 fts_bs_s_age__1  1983     -3.96 0.00494  1983
#>  5 fts_bs_s_age__1  1984     -3.90 0.00528  1984
#>  6 fts_bs_s_age__1  1985     -3.85 0.00542  1985
#>  7 fts_bs_s_age__1  1986     -3.79 0.00530  1986
#>  8 fts_bs_s_age__1  1987     -3.73 0.00509  1987
#>  9 fts_bs_s_age__1  1988     -3.67 0.00488  1988
#> 10 fts_bs_s_age__1  1989     -3.60 0.00463  1989
#> # ℹ 277 more rows

# Binomial model with cbind() response for trials data
data("qld_mortality")
# Create binomial version of the data for demonstration
qld_binomial <- qld_mortality
qld_binomial$successes <- qld_binomial$deaths
qld_binomial$failures <- qld_binomial$population - qld_binomial$deaths

mod_binomial <- ffc_gam(
  cbind(successes, failures) ~
    sex +
    fts(age,
      k = 6, bs = "cr",
      time_bs = "cr", time_k = 8
    ),
  time = "year",
  data = qld_binomial,
  family = binomial(),
  engine = "bam"
)
#> Warning: non-integer counts in a binomial glm!
summary(mod_binomial)
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> cbind(successes, failures) ~ sex + s(year, by = fts_bs_s_age__1, 
#>     bs = "cr", k = 8, m = 2, id = 1) + s(year, by = fts_bs_s_age__2, 
#>     bs = "cr", k = 8, m = 2, id = 1) + s(year, by = fts_bs_s_age__3, 
#>     bs = "cr", k = 8, m = 2, id = 1) + s(year, by = fts_bs_s_age__4, 
#>     bs = "cr", k = 8, m = 2, id = 1) + s(year, by = fts_bs_s_age__5, 
#>     bs = "cr", k = 8, m = 2, id = 1)
#> 
#> Parametric coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -5.483768   0.002294 -2390.4   <2e-16 ***
#> sexmale      0.508230   0.002166   234.6   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Approximate significance of smooth terms:
#>                           edf Ref.df Chi.sq p-value    
#> s(year):fts_bs_s_age__1 6.631  7.333  54584  <2e-16 ***
#> s(year):fts_bs_s_age__2 7.167  7.717  26231  <2e-16 ***
#> s(year):fts_bs_s_age__3 7.499  7.836  32978  <2e-16 ***
#> s(year):fts_bs_s_age__4 7.552  7.825 411279  <2e-16 ***
#> s(year):fts_bs_s_age__5 6.938  7.597 871121  <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> R-sq.(adj) =  0.982   Deviance explained = 97.7%
#> fREML =  56736  Scale est. = 1         n = 8282
```
