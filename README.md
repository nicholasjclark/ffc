
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ffc

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/ffc)](https://CRAN.R-project.org/package=ffc)
[![R-CMD-check](https://github.com/nicholasjclark/ffc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nicholasjclark/ffc/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/nicholasjclark/ffc/graph/badge.svg)](https://app.codecov.io/gh/nicholasjclark/ffc)
<!-- badges: end -->

# ffc

> **F**unctional **F**ore**C**asting

The goal of the `ffc` 📦 is to perform functional regression using
Generalized Additive Models (GAMs). The package integrates with the
extremely flexible <a href="https://cran.r-project.org/package=mgcv"
target="_blank"><code>mgcv</code></a> package to enable functional
responses to be modelled and predicted using a broad range of predictor
effects. Key among these types of predictors are *dynamic functional
predictors* using a new `fts()` term, which sets up functional
predictors whose coefficients are modelled as time-varying. These
time-varying coefficients can then be forecasted ahead using a variety
of efficient forecasting algorithms, providing unmatched flexibility to
model and predict how functional responses change over time.

## Installation

You can install the development version of ffc from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("nicholasjclark/ffc")
```

## A brief example

Load the in-built Queensland Mortality data, which contains the number
of deaths per age category over time in the state of Queensland,
Australia

``` r
library(ffc)
library(ggplot2)
theme_set(theme_bw())
data("qld_mortality")
head(qld_mortality, 15)
#>    year age    sex deaths population
#> 1  1980   0 female    190   17699.81
#> 2  1980   1 female     20   17505.27
#> 3  1980   2 female      6   17715.56
#> 4  1980   3 female      6   18080.06
#> 5  1980   4 female     10   18390.10
#> 6  1980   5 female      6   18870.54
#> 7  1980   6 female      1   19641.01
#> 8  1980   7 female      2   20475.01
#> 9  1980   8 female      2   21599.01
#> 10 1980   9 female      7   22170.09
#> 11 1980  10 female      2   21750.01
#> 12 1980  11 female      3   20866.51
#> 13 1980  12 female      1   20384.50
#> 14 1980  13 female      8   19848.04
#> 15 1980  14 female     11   19505.02
```

Visualise the observed mortality curves over time using the log10 scale

``` r
ggplot(
  data = qld_mortality,
  aes(
    x = age,
    y = deaths / population,
    group = year,
    colour = year
  )
) +
  geom_line() +
  facet_wrap(~sex) +
  scale_colour_viridis_c() +
  labs(y = "Observed log(Mortality)") +
  scale_y_log10()
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

Fit a model to estimate how the log(mortality) curve changed over time
using `deaths` as the outcome and using a time-varying function of `age`
as the primary predictor. Using `fts()`, we model the age-death function
with a set of `k = 10` thin plate basis functions whose coefficients are
allowed to vary over time, where `time = 'year'`. In this model we also
allow the time-varying effects to vary among sexes, while ensuring they
can be efficiently learned by linking their smoothing parameters. We use
the `bam()` engine (as opposed to `gam()`) for parameter estimation,
given the large size of the dataset. In future, other engines such as
`brm()` and `mvgam()`, will be made available for full luxury Bayesian
inference.

``` r
mod <- ffc_gam(
  deaths ~
    offset(log(population)) +
    sex +
    fts(
      age,
      k = 10, bs = "cr", by = sex,
      time_k = 10
    ),
  time = "year",
  data = qld_mortality,
  family = poisson(),
  engine = "bam"
)
```

Inspect the model summary; notice in the `Formula` slot how the basis
functions are modelled as `by` variables within independent smooths of
`year` that share their smoothing parameters

``` r
summary(mod)
#> 
#> Family: poisson 
#> Link function: log 
#> 
#> Formula:
#> deaths ~ sex + offset(log(population)) + s(year, by = fts_bs_s_age_bysexfemale_1, 
#>     bs = "bs", k = 10, m = 1, id = 1) + s(year, by = fts_bs_s_age_bysexfemale_2, 
#>     bs = "bs", k = 10, m = 1, id = 1) + s(year, by = fts_bs_s_age_bysexfemale_3, 
#>     bs = "bs", k = 10, m = 1, id = 1) + s(year, by = fts_bs_s_age_bysexfemale_4, 
#>     bs = "bs", k = 10, m = 1, id = 1) + s(year, by = fts_bs_s_age_bysexfemale_5, 
#>     bs = "bs", k = 10, m = 1, id = 1) + s(year, by = fts_bs_s_age_bysexfemale_6, 
#>     bs = "bs", k = 10, m = 1, id = 1) + s(year, by = fts_bs_s_age_bysexfemale_7, 
#>     bs = "bs", k = 10, m = 1, id = 1) + s(year, by = fts_bs_s_age_bysexfemale_8, 
#>     bs = "bs", k = 10, m = 1, id = 1) + s(year, by = fts_bs_s_age_bysexfemale_9, 
#>     bs = "bs", k = 10, m = 1, id = 1) + s(year, by = fts_bs_s_age_bysexmale_1, 
#>     bs = "bs", k = 10, m = 1, id = 1) + s(year, by = fts_bs_s_age_bysexmale_2, 
#>     bs = "bs", k = 10, m = 1, id = 1) + s(year, by = fts_bs_s_age_bysexmale_3, 
#>     bs = "bs", k = 10, m = 1, id = 1) + s(year, by = fts_bs_s_age_bysexmale_4, 
#>     bs = "bs", k = 10, m = 1, id = 1) + s(year, by = fts_bs_s_age_bysexmale_5, 
#>     bs = "bs", k = 10, m = 1, id = 1) + s(year, by = fts_bs_s_age_bysexmale_6, 
#>     bs = "bs", k = 10, m = 1, id = 1) + s(year, by = fts_bs_s_age_bysexmale_7, 
#>     bs = "bs", k = 10, m = 1, id = 1) + s(year, by = fts_bs_s_age_bysexmale_8, 
#>     bs = "bs", k = 10, m = 1, id = 1) + s(year, by = fts_bs_s_age_bysexmale_9, 
#>     bs = "bs", k = 10, m = 1, id = 1) + s(year, by = fts_age1_mean, 
#>     bs = "bs", k = 10, m = 1, id = 1)
#> 
#> Parametric coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -5.664337   0.004209 -1345.9   <2e-16 ***
#> sexmale      0.577793   0.005181   111.5   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Approximate significance of smooth terms:
#>                                       edf Ref.df   Chi.sq p-value    
#> s(year):fts_bs_s_age_bysexfemale_1  9.992     10  21226.4  <2e-16 ***
#> s(year):fts_bs_s_age_bysexfemale_2  9.997     10  18352.8  <2e-16 ***
#> s(year):fts_bs_s_age_bysexfemale_3  9.997     10  20256.8  <2e-16 ***
#> s(year):fts_bs_s_age_bysexfemale_4  9.999     10   8508.9  <2e-16 ***
#> s(year):fts_bs_s_age_bysexfemale_5  9.999     10    108.3  <2e-16 ***
#> s(year):fts_bs_s_age_bysexfemale_6  9.999     10  29219.5  <2e-16 ***
#> s(year):fts_bs_s_age_bysexfemale_7 10.000     10 163689.5  <2e-16 ***
#> s(year):fts_bs_s_age_bysexfemale_8 10.000     10 378156.2  <2e-16 ***
#> s(year):fts_bs_s_age_bysexfemale_9  9.999     10 318979.6  <2e-16 ***
#> s(year):fts_bs_s_age_bysexmale_1    9.995     10  34586.7  <2e-16 ***
#> s(year):fts_bs_s_age_bysexmale_2    9.999     10  23790.6  <2e-16 ***
#> s(year):fts_bs_s_age_bysexmale_3    9.999     10  26767.2  <2e-16 ***
#> s(year):fts_bs_s_age_bysexmale_4    9.999     10  11541.2  <2e-16 ***
#> s(year):fts_bs_s_age_bysexmale_5    9.999     10   1061.2  <2e-16 ***
#> s(year):fts_bs_s_age_bysexmale_6   10.000     10  58331.7  <2e-16 ***
#> s(year):fts_bs_s_age_bysexmale_7   10.000     10 231936.3  <2e-16 ***
#> s(year):fts_bs_s_age_bysexmale_8   10.000     10 394150.6  <2e-16 ***
#> s(year):fts_bs_s_age_bysexmale_9    9.996     10 125186.4  <2e-16 ***
#> s(year):fts_age1_mean               9.000      9   9489.3  <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> R-sq.(adj) =  0.986   Deviance explained = 97.7%
#> fREML =  21658  Scale est. = 1         n = 8282
```

View predicted functional curves using a fixed offset (where
`population = 1`), which allows us to calculate a standardized rate of
mortality

``` r
newdat <- qld_mortality
newdat$population <- 1
newdat$preds <- predict(
  mod,
  newdata = newdat,
  type = "response"
)

ggplot(
  data = newdat,
  aes(
    x = age,
    y = preds,
    group = year,
    colour = year
  )
) +
  geom_line() +
  facet_wrap(~sex) +
  scale_colour_viridis_c() +
  labs(y = "Expected log10(Mortality)") +
  scale_y_log10()
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

Using support from the `marginaleffects` 📦, we can make easily predict
changes in mortality rate for specific age groups. For example, here is
the expected decline in mortality rate for 17 year-olds in Queensland
over the study period

``` r
library(marginaleffects)
plot_predictions(
  mod,
  by = c("year", "sex"),
  newdata = datagrid(
    age = 17,
    year = unique,
    sex = unique,
    population = 1
  ),
  type = "response"
) +
  labs(
    x = "Year",
    y = "Expected log10(Mortality)"
  ) +
  scale_y_log10()
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

And here are the slopes of this change

``` r
plot_slopes(
  mod,
  variables = "year",
  by = c("year", "sex"),
  newdata = datagrid(
    age = 17,
    year = unique,
    sex = unique,
    population = 1
  ),
  type = "response"
) +
  labs(
    x = "Year",
    y = "1st derivative of mortality rate change"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  )
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

The time-varying coefficients can be extracted into a `tidy` format
using `fts_coefs()`, which will facilitate the use of time series models
to enable efficient forecasting of the entire curve into the future.
Using `summary = FALSE` will return draws of each coefficient time
series from the model’s empirical Bayesian posterior distribution (you
can control the number of draws that are returned using the `times`
argument):

``` r
functional_coefs <- fts_coefs(
  mod,
  summary = FALSE,
  times = 10
)
functional_coefs
#> # A tibble: 7,790 × 5
#>    .basis                     .time .estimate .realisation  year
#>    <chr>                      <int>     <dbl>        <int> <int>
#>  1 fts_bs_s_age_bysexfemale_1  1980     -3.78            1  1980
#>  2 fts_bs_s_age_bysexfemale_1  1981     -3.79            1  1981
#>  3 fts_bs_s_age_bysexfemale_1  1982     -3.80            1  1982
#>  4 fts_bs_s_age_bysexfemale_1  1983     -3.82            1  1983
#>  5 fts_bs_s_age_bysexfemale_1  1984     -3.83            1  1984
#>  6 fts_bs_s_age_bysexfemale_1  1985     -3.83            1  1985
#>  7 fts_bs_s_age_bysexfemale_1  1986     -3.81            1  1986
#>  8 fts_bs_s_age_bysexfemale_1  1987     -3.80            1  1987
#>  9 fts_bs_s_age_bysexfemale_1  1988     -3.78            1  1988
#> 10 fts_bs_s_age_bysexfemale_1  1989     -3.76            1  1989
#> # ℹ 7,780 more rows
```

The basis function coefficient time series can be plotted using
`autoplot()`

``` r
autoplot(functional_coefs)
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

Clearly there is a lot of structure and dependence here, suggesting that
a dynamic factor model fitted to these coefficient time series would be
valuable for creating functional forecasts. But for now we can apply any
model from the `fable` 📦 to these replicate time series and generate
future forecast realisations, which can be summarised to approximate the
full uncertainty in our coefficient forecast distributions. Again here
you can control the number of forecast paths that are simulated from the
underlying time series models using the `times` argument

``` r
functional_fc <- forecast(
  object = functional_coefs,
  h = 5,
  times = 5
)
#> Registered S3 method overwritten by 'tsibble':
#>   method               from 
#>   as_tibble.grouped_df dplyr
functional_fc
#> # A tsibble: 4,750 x 6 [1Y]
#> # Key:       .basis, .realisation, .model, .rep [950]
#>    .basis        .realisation .model .time .rep    .sim
#>    <chr>                <int> <chr>  <dbl> <chr>  <dbl>
#>  1 fts_age1_mean            1 ARIMA   2021 1     -0.384
#>  2 fts_age1_mean            1 ARIMA   2022 1     -0.393
#>  3 fts_age1_mean            1 ARIMA   2023 1     -0.403
#>  4 fts_age1_mean            1 ARIMA   2024 1     -0.415
#>  5 fts_age1_mean            1 ARIMA   2025 1     -0.434
#>  6 fts_age1_mean            1 ARIMA   2021 2     -0.387
#>  7 fts_age1_mean            1 ARIMA   2022 2     -0.405
#>  8 fts_age1_mean            1 ARIMA   2023 2     -0.428
#>  9 fts_age1_mean            1 ARIMA   2024 2     -0.448
#> 10 fts_age1_mean            1 ARIMA   2025 2     -0.467
#> # ℹ 4,740 more rows
```

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/nicholasjclark/ffc/issues)

## Contributing

Contributions are very welcome, but please see our [Code of
Conduct](https://github.com/nicholasjclark/ffc/blob/main/.github/CODE_OF_CONDUCT.md)
when you are considering changes that you would like to make.

## License

The `ffc` project is licensed under an `MIT` open source license
