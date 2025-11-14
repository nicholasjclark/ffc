# Modeling Mortality Trends with Functional Forecasting

``` r
library(ffc)
library(ggplot2)
library(marginaleffects)
theme_set(theme_bw())
```

## Introduction

Mortality patterns represent one of the most fundamental demographic
processes, yet traditional statistical approaches often miss the
complex, evolving functional relationships between age and death rates
over time. This vignette demonstrates how functional forecasting using
the `ffc` package captures these dynamic relationships that change
shape—not just magnitude—over time.

### What makes mortality data ideal for functional forecasting?

1.  Complex age patterns: Mortality exhibits characteristic J-shaped
    curves across age groups
2.  Temporal evolution: These functional relationships shift
    systematically over time  
3.  Smooth transitions: Changes occur gradually, making them suitable
    for GAM-based smoothing
4.  Forecasting relevance: Understanding future mortality trends has
    critical policy implications

### Key concepts we’ll explore

- Time-varying coefficients: How functional relationships evolve over
  time
- Hierarchical smoothing: Modeling shared trends with group-specific
  deviations  
- Functional forecasting: Predicting entire curves into the future
- Model diagnostics: Validating functional time series models

## Data exploration

We’ll use Queensland mortality data spanning 1980-2020, containing death
counts by age, sex, and year with corresponding population denominators.

``` r
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

The dataset structure is ideal for functional analysis: we have a
functional predictor (age) whose relationship with the outcome
(mortality) changes over time (year), with hierarchical structure (sex).

### Visualizing mortality patterns

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
  geom_line(alpha = 0.7) +
  facet_wrap(~sex) +
  scale_colour_viridis_c(name = "Year") +
  labs(
    x = "Age", 
    y = "Mortality rate",
    title = "Queensland mortality patterns over time"
  ) +
  scale_y_log10()
```

![Observed mortality rates by age in Queensland, 1980-2020. The
characteristic J-shaped curves show systematic downward shifts over
time, indicating mortality improvements across all age
groups.](mortality-analysis_files/figure-html/mortality-curves-1.png)

Observed mortality rates by age in Queensland, 1980-2020. The
characteristic J-shaped curves show systematic downward shifts over
time, indicating mortality improvements across all age groups.

Key observations:

- J-shaped curves: High infant mortality, low childhood mortality,
  exponential increase with age
- Temporal shifts: Consistent downward movement of entire curves over
  time  
- Sex differences: Males consistently higher mortality, similar temporal
  patterns
- Functional evolution: The shape itself evolves, not just vertical
  shifts
- Overall decline: Mortality rates have declined substantially across
  all ages, reflecting major improvements in healthcare, lifestyle, and
  living conditions

## Functional modeling approach

### Why traditional models fall short

Standard approaches might model this as:

``` r
# Traditional approach - misses functional evolution
glm(deaths ~ age + year + sex, family = poisson(), offset = log(population))
```

This assumes linear age effects and additive time trends—clearly
inadequate for evolving functional relationships.

### The ffc solution: time-varying coefficients

The [`fts()`](https://nicholasjclark.github.io/ffc/reference/fts.md)
function creates basis functions whose coefficients evolve over time,
capturing how the age-mortality relationship changes:

``` r
mod <- ffc_gam(
  deaths ~
    offset(log(population)) +
    sex +
    # Time-varying level: shared temporal trends
    fts(
      year,
      mean_only = TRUE,
      bs = "tp",
      time_k = 35,  # Temporal resolution
      time_m = 1    # Smoothness constraint
    ) +
    # Time-varying age effects: sex-specific deviations
    fts(
      age,
      by = sex,
      bs = "tp", 
      time_k = 15,  # Lower temporal resolution for age effects
      time_m = 1
    ),
  time = "year",
  data = qld_mortality,
  family = poisson(),
  engine = "bam"  # Efficient for large datasets
)
```

### Understanding the model structure

``` r
summary(mod)
#> 
#> Family: poisson 
#> Link function: log 
#> 
#> Formula:
#> deaths ~ sex + offset(log(population)) + s(year, by = fts_year1_mean, 
#>     bs = "ts", k = 35, m = 1, id = 1) + s(year, by = fts_bs_s_age_bysexfemale_1, 
#>     bs = "ts", k = 15, m = 1, id = 2) + s(year, by = fts_bs_s_age_bysexfemale_2, 
#>     bs = "ts", k = 15, m = 1, id = 2) + s(year, by = fts_bs_s_age_bysexfemale_3, 
#>     bs = "ts", k = 15, m = 1, id = 2) + s(year, by = fts_bs_s_age_bysexfemale_4, 
#>     bs = "ts", k = 15, m = 1, id = 2) + s(year, by = fts_bs_s_age_bysexfemale_5, 
#>     bs = "ts", k = 15, m = 1, id = 2) + s(year, by = fts_bs_s_age_bysexfemale_6, 
#>     bs = "ts", k = 15, m = 1, id = 2) + s(year, by = fts_bs_s_age_bysexfemale_7, 
#>     bs = "ts", k = 15, m = 1, id = 2) + s(year, by = fts_bs_s_age_bysexfemale_8, 
#>     bs = "ts", k = 15, m = 1, id = 2) + s(year, by = fts_bs_s_age_bysexfemale_9, 
#>     bs = "ts", k = 15, m = 1, id = 2) + s(year, by = fts_bs_s_age_bysexmale_1, 
#>     bs = "ts", k = 15, m = 1, id = 2) + s(year, by = fts_bs_s_age_bysexmale_2, 
#>     bs = "ts", k = 15, m = 1, id = 2) + s(year, by = fts_bs_s_age_bysexmale_3, 
#>     bs = "ts", k = 15, m = 1, id = 2) + s(year, by = fts_bs_s_age_bysexmale_4, 
#>     bs = "ts", k = 15, m = 1, id = 2) + s(year, by = fts_bs_s_age_bysexmale_5, 
#>     bs = "ts", k = 15, m = 1, id = 2) + s(year, by = fts_bs_s_age_bysexmale_6, 
#>     bs = "ts", k = 15, m = 1, id = 2) + s(year, by = fts_bs_s_age_bysexmale_7, 
#>     bs = "ts", k = 15, m = 1, id = 2) + s(year, by = fts_bs_s_age_bysexmale_8, 
#>     bs = "ts", k = 15, m = 1, id = 2) + s(year, by = fts_bs_s_age_bysexmale_9, 
#>     bs = "ts", k = 15, m = 1, id = 2)
#> 
#> Parametric coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -5.648061   0.004039 -1398.2   <2e-16 ***
#> sexmale      0.574335   0.004986   115.2   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Approximate significance of smooth terms:
#>                                      edf Ref.df Chi.sq p-value    
#> s(year):fts_year1_mean             32.05     34  10342  <2e-16 ***
#> s(year):fts_bs_s_age_bysexfemale_1 14.50     15  16349  <2e-16 ***
#> s(year):fts_bs_s_age_bysexfemale_2 12.85     15  10267  <2e-16 ***
#> s(year):fts_bs_s_age_bysexfemale_3 14.49     15   9610  <2e-16 ***
#> s(year):fts_bs_s_age_bysexfemale_4 13.48     15   8470  <2e-16 ***
#> s(year):fts_bs_s_age_bysexfemale_5 14.41     15   8433  <2e-16 ***
#> s(year):fts_bs_s_age_bysexfemale_6 13.35     15   8667  <2e-16 ***
#> s(year):fts_bs_s_age_bysexfemale_7 14.35     15   7764  <2e-16 ***
#> s(year):fts_bs_s_age_bysexfemale_8 10.12     15  11476  <2e-16 ***
#> s(year):fts_bs_s_age_bysexfemale_9 13.05     15   6036  <2e-16 ***
#> s(year):fts_bs_s_age_bysexmale_1   14.44     15  19549  <2e-16 ***
#> s(year):fts_bs_s_age_bysexmale_2   12.95     15  11146  <2e-16 ***
#> s(year):fts_bs_s_age_bysexmale_3   14.64     15  12750  <2e-16 ***
#> s(year):fts_bs_s_age_bysexmale_4   13.75     15  11841  <2e-16 ***
#> s(year):fts_bs_s_age_bysexmale_5   14.56     15  14649  <2e-16 ***
#> s(year):fts_bs_s_age_bysexmale_6   13.74     15  12002  <2e-16 ***
#> s(year):fts_bs_s_age_bysexmale_7   14.52     15  11843  <2e-16 ***
#> s(year):fts_bs_s_age_bysexmale_8   10.19     15  13491  <2e-16 ***
#> s(year):fts_bs_s_age_bysexmale_9   12.44     15   6716  <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> R-sq.(adj) =  0.986   Deviance explained = 97.6%
#> fREML =  22189  Scale est. = 1         n = 8282
```

Key components:

1.  Fixed effects: `sex` captures baseline male-female differences
2.  Time-varying level: `fts(year, mean_only=TRUE)` models shared
    temporal trends  
3.  Time-varying functions: `fts(age, by=sex)` captures how age patterns
    evolve differently by sex

The model automatically creates a hierarchical structure where basis
functions become `by` variables in smooths of time, sharing smoothing
parameters for computational efficiency.

## Model interpretation

### Predicted functional curves

``` r
newdat <- qld_mortality
newdat$population <- 1  # Standardized predictions

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
  scale_colour_viridis_c(name = "Year") +
  labs(
    x = "Age",
    y = "Expected mortality rate",
    title = "Model predictions: evolving mortality patterns"
  ) +
  scale_y_log10()
```

![Model predictions showing expected mortality curves by age and year.
Smooth curves demonstrate systematic evolution of the age-mortality
relationship over four
decades.](mortality-analysis_files/figure-html/predicted-curves-1.png)

Model predictions showing expected mortality curves by age and year.
Smooth curves demonstrate systematic evolution of the age-mortality
relationship over four decades.

The model successfully captures the smooth evolution of J-shaped
mortality curves while preserving the characteristic functional form.
The downward trend across all mortality curves demonstrates substantial
improvement in survival rates over the 40-year period. This reflects
major advances in healthcare, improved living conditions, better
nutrition, and reduced exposure to environmental hazards. Notably, the
decline is not uniform across all ages - the model captures how the rate
of improvement varies by age group and sex, with some periods showing
more rapid mortality decline than others.

### Focused analysis: teenage mortality trends

Let’s examine how mortality patterns changed for a specific age group:

``` r
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
    y = "Expected mortality rate",
    title = "Teenage mortality trends over time"
  ) +
  scale_y_log10()
```

![Predicted mortality trends for 17-year-olds in Queensland, 1980-2020.
Shows approximately 70% decline in mortality rates with consistent
male-female
differences.](mortality-analysis_files/figure-html/mortality-trends-1.png)

Predicted mortality trends for 17-year-olds in Queensland, 1980-2020.
Shows approximately 70% decline in mortality rates with consistent
male-female differences.

### Rate of change analysis

Understanding *how fast* mortality is improving:

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
    y = "Rate of mortality change",
    title = "Speed of mortality improvement over time"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    alpha = 0.7
  )
```

![First derivatives showing the rate of mortality improvement for
17-year-olds. Fluctuations reveal periods of faster and slower mortality
decline.](mortality-analysis_files/figure-html/mortality-derivatives-1.png)

First derivatives showing the rate of mortality improvement for
17-year-olds. Fluctuations reveal periods of faster and slower mortality
decline.

## Model diagnostics

### Residual analysis

``` r
# Extract residuals
resids <- residuals(mod)
fitted_vals <- fitted(mod)

# Create data frame for ggplot
resid_data <- data.frame(
  fitted = fitted_vals,
  residuals = resids
)

# ggplot version
ggplot(resid_data, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Fitted values", 
    y = "Residuals",
    title = "Residual analysis"
  ) +
  theme_bw()
```

![Model residuals showing good fit with no systematic patterns. Random
scatter around zero indicates appropriate model
specification.](mortality-analysis_files/figure-html/residuals-1.png)

Model residuals showing good fit with no systematic patterns. Random
scatter around zero indicates appropriate model specification.

### Model diagnostics

For functional time series models, traditional GAM diagnostics need to
be adapted. The ffc package will develop specialized diagnostic methods
for:

- **Basis adequacy**: Ensuring sufficient flexibility for time-varying
  coefficients
- **Temporal autocorrelation**: Checking coefficient time series for
  proper structure  
- **Functional residuals**: Examining residuals at both observation and
  functional levels
- **Forecast validation**: Cross-validation approaches for functional
  forecasting

Development of comprehensive diagnostic tools for ffc models is ongoing.

## Extracting time-varying coefficients

The power of functional forecasting lies in treating the evolving
coefficients as time series that can be forecasted:

``` r
functional_coefs <- fts_coefs(
  mod,
  summary = FALSE,
  times = 10  # Number of coefficient draws
)
functional_coefs
#> # A tibble: 7,790 × 5
#>    .basis         .time .estimate .realisation  year
#>  * <chr>          <int>     <dbl>        <int> <int>
#>  1 fts_year1_mean  1980     0.387            1  1980
#>  2 fts_year1_mean  1981     0.377            1  1981
#>  3 fts_year1_mean  1982     0.415            1  1982
#>  4 fts_year1_mean  1983     0.319            1  1983
#>  5 fts_year1_mean  1984     0.309            1  1984
#>  6 fts_year1_mean  1985     0.330            1  1985
#>  7 fts_year1_mean  1986     0.280            1  1986
#>  8 fts_year1_mean  1987     0.253            1  1987
#>  9 fts_year1_mean  1988     0.240            1  1988
#> 10 fts_year1_mean  1989     0.272            1  1989
#> # ℹ 7,780 more rows
```

### Visualizing coefficient evolution

``` r
autoplot(functional_coefs) +
  labs(title = "Evolution of functional basis coefficients")
```

![Time series of functional basis coefficients. Each panel shows how
different aspects of the age-mortality relationship evolved over time,
revealing complex temporal
dependencies.](mortality-analysis_files/figure-html/plot-coefficients-1.png)

Time series of functional basis coefficients. Each panel shows how
different aspects of the age-mortality relationship evolved over time,
revealing complex temporal dependencies.

Interpretation:

- Complex patterns: Each coefficient series shows distinct temporal
  behavior
- Smooth evolution: Changes occur gradually, suitable for time series
  modeling
- Interdependence: Multiple coefficients work together to create the
  functional evolution

## Forecasting mortality patterns

### Coefficient forecasting

``` r
functional_fc <- forecast(
  object = functional_coefs,
  h = 5,      # 5-year horizon
  times = 10, # Forecast replicates
  model = "ARIMA"
)
functional_fc
#> # A tsibble: 9,500 x 6 [1Y]
#> # Key:       .basis, .realisation, .model, .rep [1,900]
#>    .basis                     .realisation .model .time .rep   .sim
#>    <chr>                             <int> <chr>  <dbl> <chr> <dbl>
#>  1 fts_bs_s_age_bysexfemale_1            1 ARIMA   2021 1      4.68
#>  2 fts_bs_s_age_bysexfemale_1            1 ARIMA   2022 1      4.70
#>  3 fts_bs_s_age_bysexfemale_1            1 ARIMA   2023 1      4.75
#>  4 fts_bs_s_age_bysexfemale_1            1 ARIMA   2024 1      4.84
#>  5 fts_bs_s_age_bysexfemale_1            1 ARIMA   2025 1      4.97
#>  6 fts_bs_s_age_bysexfemale_1            1 ARIMA   2021 10     4.69
#>  7 fts_bs_s_age_bysexfemale_1            1 ARIMA   2022 10     4.75
#>  8 fts_bs_s_age_bysexfemale_1            1 ARIMA   2023 10     4.85
#>  9 fts_bs_s_age_bysexfemale_1            1 ARIMA   2024 10     4.97
#> 10 fts_bs_s_age_bysexfemale_1            1 ARIMA   2025 10     5.05
#> # ℹ 9,490 more rows
```

The forecast includes uncertainty in both the time series models and the
original coefficient estimation, providing realistic prediction
intervals.

### Future mortality patterns

Let’s demonstrate how to generate complete mortality curve forecasts by
combining the forecasted coefficients with the model structure:

``` r
# Create forecast data for future years
future_years <- 2021:2025
newdata_forecast <- expand.grid(
  age = unique(qld_mortality$age),
  sex = unique(qld_mortality$sex),
  year = future_years,
  population = 1
)

# Generate forecasted mortality curves using forecast()
# This integrates the time series forecasts of coefficients
mortality_forecasts <- forecast(
  object = mod,
  newdata = newdata_forecast,
  model = "ARIMA",
  type = "expected"
)
head(mortality_forecasts)
#> # A tibble: 6 × 6
#>   .estimate  .error    .q2.5     .q10     .q90   .q97.5
#>       <dbl>   <dbl>    <dbl>    <dbl>    <dbl>    <dbl>
#> 1  0.00178  0.00166 0.00152  0.00162  0.00193  0.00203 
#> 2  0.00113  0.00111 0.000987 0.00104  0.00122  0.00127 
#> 3  0.000726 0.00131 0.000624 0.000668 0.000775 0.000807
#> 4  0.000468 0.00149 0.000409 0.000431 0.000498 0.000524
#> 5  0.000306 0.00165 0.000266 0.000284 0.000327 0.000346
#> 6  0.000206 0.00186 0.000178 0.000190 0.000223 0.000231

# Plot forecasted mortality rates, together with uncertainties
ggplot(mortality_forecasts %>%
         dplyr::bind_cols(newdata_forecast),
       aes(x = age,
           y = .estimate,
           group = year,
           colour = year)) +
  # geom_ribbon(aes(
  #   ymin = .q10.0,
  #   ymax = .q90.0,
  #   fill = year
  # ),
  # alpha = 0.2,
  # colour = NA) +
  geom_line() +
  facet_wrap(~sex) +
  scale_colour_viridis_c(name = "Year") +
  # scale_fill_viridis_c(name = "Year") +
  labs(
    x = "Age", 
    y = "Mortality rate",
    title = "Expected mortality"
  ) +
  scale_y_log10()
```

![Forecasted mortality curves for Queensland, 2021-2025. Shows how the
functional forecasting approach predicts future age-mortality
patterns.](mortality-analysis_files/figure-html/forecast-curves-1.png)

Forecasted mortality curves for Queensland, 2021-2025. Shows how the
functional forecasting approach predicts future age-mortality patterns.

This forecasting approach enables prediction of entire functional
relationships by combining time series forecasts of the functional
coefficients with the underlying model structure.

## Key insights and methodology

### What we learned about Queensland mortality

- Systematic improvement: Mortality declined across all age groups over
  40 years
- Functional evolution: Changes involved curve shape, not just magnitude
  shifts
- Heterogeneous trends: Rate of improvement varied by age and sex
- Complex patterns: Simple parametric models would miss these
  relationships

### Why functional forecasting matters

1.  Captures complexity: Models evolving functional relationships
    traditional methods miss
2.  Preserves structure: Maintains biological/physical constraints in
    forecasts  
3.  Quantifies uncertainty: Propagates uncertainty through the entire
    functional evolution
4.  Policy relevance: Enables sophisticated demographic projections

### Model design principles

- Hierarchical structure: Shared trends with group-specific deviations
- Temporal smoothing: Balance flexibility with smoothness using `time_k`
  and `time_m`
- Basis choice: Thin plate splines work well for age and time
- Computational efficiency:
  [`bam()`](https://nicholasjclark.github.io/ffc/reference/bam.md)
  engine handles large datasets effectively

## Extensions and further analysis

### Additional applications

This methodology extends to many functional forecasting problems: -
**Environmental data**: Temperature/precipitation profiles over time -
**Economic indicators**: Yield curves in finance  
- **Biological processes**: Growth curves, dose-response relationships -
**Spatial processes**: Evolving spatial fields

### Advanced features

The `ffc` package supports additional complexity: - **Multiple
functional predictors**: `fts(age) + fts(time_since_event)`  
- **Interaction effects**: Functional relationships that depend on other
variables - **Different time scales**: Multiple temporal resolutions in
the same model - **Alternative distributions**: Beyond Poisson to handle
different data types

## Conclusion

Functional forecasting with `ffc` provides a principled approach to
modeling and predicting evolving functional relationships. By treating
time-varying coefficients as forecastable time series, we can:

1.  **Capture complexity** that traditional methods miss
2.  **Generate realistic forecasts** of entire functional
    relationships  
3.  **Quantify uncertainty** appropriately across all sources
4.  **Provide actionable insights** for policy and planning

The Queensland mortality analysis demonstrates these capabilities in a
real-world context, revealing patterns and trends that inform our
understanding of demographic change and enable sophisticated projection
methods for public health planning.
