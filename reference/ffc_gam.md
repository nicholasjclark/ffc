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
#> 3.985387 3.572323 3.185908 2.815229 2.452154 2.112488 
head(predict(mod, type = "response"))
#>         1         2         3         4         5         6 
#> 53.806134 35.599209 24.189249 16.697007 11.613336  8.268789 

# Extract basis coefficient time series
(functional_ts <- fts_coefs(mod))
#> # A tibble: 287 × 5
#>    .basis          .time .estimate     .se  year
#>  * <chr>           <int>     <dbl>   <dbl> <int>
#>  1 fts_bs_s_age__1  1980     -4.12 0.00827  1980
#>  2 fts_bs_s_age__1  1981     -4.07 0.00655  1981
#>  3 fts_bs_s_age__1  1982     -4.01 0.00538  1982
#>  4 fts_bs_s_age__1  1983     -3.96 0.00490  1983
#>  5 fts_bs_s_age__1  1984     -3.90 0.00489  1984
#>  6 fts_bs_s_age__1  1985     -3.84 0.00497  1985
#>  7 fts_bs_s_age__1  1986     -3.79 0.00498  1986
#>  8 fts_bs_s_age__1  1987     -3.73 0.00501  1987
#>  9 fts_bs_s_age__1  1988     -3.66 0.00511  1988
#> 10 fts_bs_s_age__1  1989     -3.60 0.00519  1989
#> # ℹ 277 more rows

# Binomial model with cbind() response for trials data
data("qld_mortality")
# Create binomial version of the data for demonstration
qld_binomial <- qld_mortality
qld_binomial$trials <- qld_binomial$population
qld_binomial$successes <- qld_binomial$deaths

mod_binomial <- ffc_gam(
  cbind(successes, trials - successes) ~
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
#> Error: {var_type_cap} {.field {missing_vars}} not found in data
summary(mod_binomial)
#> Error: object 'mod_binomial' not found
```
