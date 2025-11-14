# Convert ffc_gam forecasts to fable object

Convert ffc_gam forecasts to fable object

## Usage

``` r
# S3 method for class 'ffc_gam'
as_fable(
  object,
  newdata,
  forecasts = NULL,
  response = NULL,
  model = "ARIMA",
  key_vars = NULL,
  ...
)
```

## Arguments

- object:

  An `ffc_gam` object that has been used to generate forecasts

- newdata:

  A `data.frame` containing the forecast data with the response variable

- forecasts:

  Optional pre-computed forecasts from
  [`forecast.ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/forecast.ffc_gam.md).
  If NULL, forecasts will be generated using the remaining arguments

- response:

  Character string specifying the response variable name. If NULL,
  automatically detected from the model formula

- model:

  A character string representing the forecasting model to use if
  generating forecasts. Default is "ARIMA"

- key_vars:

  Optional character vector specifying grouping variables. If NULL,
  automatically detected from categorical variables in newdata

- ...:

  Additional arguments passed to
  [`forecast.ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/forecast.ffc_gam.md)
  if generating forecasts

## Value

A `fable` object containing forecast distributions and point estimates

## Details

Converts forecasting output from `ffc_gam` objects into a properly
formatted `fable` object that can be used with `fabletools` functions
like
[`autoplot()`](https://nicholasjclark.github.io/ffc/reference/autoplot.md),
[`accuracy()`](https://generics.r-lib.org/reference/accuracy.html), and
forecast combination methods

## See also

[`forecast.ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/forecast.ffc_gam.md)

## Author

Nicholas J Clark

## Examples

``` r
# \donttest{
# Basic usage with automatic detection
library(fable)
#> Loading required package: fabletools
library(tsibble)
#> 
#> Attaching package: ‘tsibble’
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, union
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union

# Prepare tourism data
tourism_melb <- tourism %>%
  filter(Region == "Melbourne", Purpose == "Visiting") %>%
  mutate(quarter = as.numeric(substr(as.character(Quarter), 6, 6)), 
         time = row_number())
#> Warning: There was 1 warning in `mutate()`.
#> ℹ In argument: `quarter = as.numeric(substr(as.character(Quarter), 6, 6))`.
#> Caused by warning:
#> ! NAs introduced by coercion

# Split data
train <- tourism_melb %>% slice_head(n = 75)
test <- tourism_melb %>% slice_tail(n = 5)

# Fit model
mod <- ffc_gam(
  Trips ~ fts(time, mean_only = TRUE, time_k = 50, time_m = 1) +
          fts(quarter, k = 4, time_k = 15, time_m = 1),
  time = "time", data = train, family = tw(), engine = "gam"
)
#> Error in init_gam(formula(formula), data = dat, family = family): Not enough (non-NA) data to do anything meaningful

# Convert to fable with auto-detection
fc_fable <- as_fable(mod, newdata = test, model = "ETS")
#> Error: object 'mod' not found

# Use fabletools ecosystem
autoplot(fc_fable, train)  # Forecast plot
#> Error: object 'fc_fable' not found
accuracy(fc_fable, test)   # Accuracy metrics
#> Error: object 'fc_fable' not found
hilo(fc_fable, level = c(80, 95))  # Prediction intervals
#> Error: object 'fc_fable' not found

# Distribution summaries
fc_fable %>%
  summarise(
    mean_forecast = mean(.dist),
    q25 = quantile(.dist, 0.25),
    q75 = quantile(.dist, 0.75)
  )
#> Error: object 'fc_fable' not found

# With pre-computed forecasts
forecasts <- forecast(mod, newdata = test, summary = FALSE)
#> Error: object 'mod' not found
fc_fable2 <- as_fable(mod, newdata = test, forecasts = forecasts)
#> Error: object 'mod' not found

# With custom parameters
fc_fable3 <- as_fable(
  mod, 
  newdata = test, 
  model = "ARIMA",
  response = "Trips",
  key_vars = c("Region", "State")
)
#> Error: object 'mod' not found

# Model comparison workflow
fc_arima <- as_fable(mod, newdata = test, model = "ARIMA")
#> Error: object 'mod' not found
fc_ets <- as_fable(mod, newdata = test, model = "ETS")
#> Error: object 'mod' not found

# Combine and compare
combined <- bind_rows(fc_arima, fc_ets)
#> Error: object 'fc_arima' not found
autoplot(combined, train)
#> Error: object 'combined' not found
accuracy(combined, test)
#> Error: object 'combined' not found
# }
```
