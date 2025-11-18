# Fit a Vector autoregressive dynamic factor model

Fit a Vector autoregressive dynamic factor model using Stan

## Usage

``` r
VARDF(
  formula,
  family = gaussian(),
  h = get_stan_param("h", "forecast"),
  chains = get_stan_param("chains"),
  cores = get_stan_param("cores"),
  iter = get_stan_param("iter"),
  warmup = floor(iter/2),
  adapt_delta = get_stan_param("adapt_delta"),
  max_treedepth = get_stan_param("max_treedepth"),
  silent = get_stan_param("silent"),
  ...
)
```

## Arguments

- formula:

  Model specification (see "Specials" section)

- family:

  A family object specifying the outcome distribution to use in fitting.
  Currently only [`gaussian()`](https://rdrr.io/r/stats/family.html) and
  [`scat()`](https://rdrr.io/pkg/mgcv/man/scat.html) (i.e. Student-T)
  are supported

- h:

  `integer` specifying the forecast horizon

- chains:

  `integer` specifying the number of chains to be run

- cores:

  `integer` specifying the number of parallel cores to use

- iter:

  `integer` specifying the total number of iterations to run per chain
  (including warmup)

- warmup:

  `integer` specifying the number of initial iterations to use as burnin

- adapt_delta:

  the thin of the jumps in a HMC method

- max_treedepth:

  maximum tree depth per iteration

- silent:

  Logical indicating whether to suppress Stan sampling progress output.
  Default is `TRUE`. When `FALSE`, shows progress approximately every
  10% of iterations.

- ...:

  other arguments to pass to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)

## Value

A model specification

## Author

Nicholas J Clark

## Examples

``` r
# Fit a functional forecasting model, then use VARDF for forecasting
library(dplyr)

# Split growth data into training and test sets
train_data <- growth_data |> filter(age_yr <= 13)
test_data <- growth_data |> filter(age_yr > 13)

# Step 1: Fit ffc_gam model with time-varying coefficients
mod <- ffc_gam(
  height_cm ~
    s(id, bs = "re") +
    fts(age_yr, time_k = 5),
  data = train_data,
  time = "age_yr",
  family = gaussian()
)

# Step 2: Use VARDF for forecasting functional coefficients
fc <- forecast(mod, newdata = test_data, model = "VARDF",
               chains = 1, iter = 300)
```
