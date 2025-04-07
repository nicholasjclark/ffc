# Working example to show how ffc_gam() can be used
library(vital)
library(tidyverse); theme_set(theme_classic(base_size = 12))
library(mgcv)
library(marginaleffects)

# Grab some functional data from the vital package
dat <- aus_mortality %>%
  dplyr::filter(Year > 1970,
                State == 'Queensland',
                Sex != 'total') %>%
  dplyr::mutate(Log_population = log(Exposure),
                Zeros = ifelse(Age == 0, 1, 0),
                Deaths = ceiling(Deaths),
                Sex = factor(Sex))

# Fit a model to estimate how rates of mortality change over time
mod <- ffc_gam(
  formula = Deaths ~
    offset(log_population) +
    Sex * Zeros +

    # Use 6 basis functions to characterise the age-mortality
    # function; for each of these functions, use a thin-plate
    # spline to estimate how it's coefficient changes through time
    fts(Age, bs = 'tp', k = 6,
        time_k = 25, time_bs = 'tp', by = Sex),
  data = dat,
  time = 'Year',
  family = poisson(),
  engine = 'bam',
  discrete = TRUE
)

# Have a look at the model object
summary(mod)
class(mod)
mod$call
mod$fts_smooths
mod$gam_init

# Inspect hindcast predictions of log(mortality)
plot_predictions(
  mod,
  by = c('Age', 'Year', 'Sex'),
  type = 'link'
) +
  labs(y = 'log(mortality) rate',
       title = 'Standardised mortality rates')

dat$preds <- predict(mod, newdata = dat, type = 'link')

ggplot(dat,
       aes(x = Age,
           y = preds,
           colour = Year,
           group = Year)) +
  geom_smooth(se = FALSE) +
  facet_wrap(~ Sex) +
  scale_colour_viridis_c() +
  labs(y = 'log(mortality) rate',
       title = 'Standardised mortality rates')

# To do:
# 1. Add functionality to extract the basis function coefficient time series

# Extract basis coefficients for forecasting by using
# predict(..., type = 'terms')
predict(mod,
        newdata = dat,
        type = 'terms')

# 2. Add functionality to forecast the basis function coefficients, and to
#    combine these with the rest of the betas so that forecasts can be computed
# 3. Write unit tests and ensure the functions work with a broad range of bases,
#   including MRF smooths, cyclic smooths etc...
