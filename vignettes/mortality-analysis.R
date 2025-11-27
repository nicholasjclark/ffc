## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)


## ----setup--------------------------------------------------------------------
library(ffc)
library(ggplot2)
library(marginaleffects)
theme_set(theme_bw())


## -----------------------------------------------------------------------------
data("qld_mortality")
head(qld_mortality, 15)


## ----mortality-curves, fig.cap="Observed mortality rates by age in Queensland, 1980-2020. The characteristic J-shaped curves show systematic downward shifts over time, indicating mortality improvements across all age groups."----
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


## ----fit-model, fig.show='hide'-----------------------------------------------
mod <- ffc_gam(
  deaths ~
    offset(log(population)) +
    sex +
    age +
    # Time-varying level: shared temporal trends
    fts(
      year,
      mean_only = TRUE,
      bs = "tp",
      time_k = 35,
      time_m = 1
    ) +
    # Time-varying age effects: sex-specific deviations
    fts(
      age,
      by = sex,
      bs = "tp",
      time_k = 15,
      time_m = 1
    ),
  time = "year",
  data = qld_mortality,
  family = poisson(),
  # Efficient for large datasets
  engine = "bam"
)


## ----model-summary------------------------------------------------------------
summary(mod)


## ----predicted-curves, fig.cap="Model predictions showing expected mortality curves by age and year. Smooth curves demonstrate systematic evolution of the age-mortality relationship over four decades."----
newdat <- qld_mortality
newdat$population <- 1 # Standardized predictions

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


## ----mortality-trends, fig.cap="Predicted mortality trends for 17-year-olds in Queensland, 1980-2020. Shows approximately 70% decline in mortality rates with consistent male-female differences."----
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


## ----mortality-derivatives, fig.cap="First derivatives showing the rate of mortality improvement for 17-year-olds. Fluctuations reveal periods of faster and slower mortality decline."----
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


## ----extract-coefficients-----------------------------------------------------
functional_coefs <- fts_coefs(
  mod,
  summary = FALSE,
  times = 10
)
functional_coefs


## ----plot-coefficients, fig.cap="Time series of functional basis coefficients. Each panel shows how different aspects of the age-mortality relationship evolved over time, revealing complex temporal dependencies."----
autoplot(functional_coefs) +
  labs(title = "Evolution of functional basis coefficients")


## ----forecast-coefficients----------------------------------------------------
functional_fc <- forecast(
  object = functional_coefs,
  h = 5, # 5-year horizon
  times = 10, # Forecast replicates
  model = "ARIMA"
)
functional_fc


## ----forecast-curves, fig.cap="Forecasted mortality curves for Queensland, 2021-2025. Shows how the functional forecasting approach predicts future age-mortality patterns."----
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
  mean_model = "ENS",
  type = "expected"
)
head(mortality_forecasts)

# Plot forecasted mortality rates, together with uncertainties
ggplot(
  mortality_forecasts |>
    dplyr::bind_cols(newdata_forecast),
  aes(x = age, y = .estimate, group = year, colour = year)
) +
  geom_ribbon(
    aes(
      ymin = .q10,
      ymax = .q90,
      fill = year
    ),
    alpha = 0.2,
    colour = NA
  ) +
  geom_line() +
  facet_wrap(~sex) +
  scale_colour_viridis_c(name = "Year") +
  scale_fill_viridis_c(name = "Year") +
  labs(
    x = "Age",
    y = "Mortality rate",
    title = "Expected mortality"
  ) +
  scale_y_log10()
