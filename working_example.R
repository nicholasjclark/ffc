# Working example to show how ffc_gam() can be used
library(vital)
library(tidyverse); theme_set(theme_bw(base_size = 12))
library(mgcv)
library(gratia)
library(marginaleffects)

# Grab some functional data from the vital package
dat <- aus_mortality %>%
  dplyr::filter(Year > 1950,
                State == 'Queensland',
                Sex != 'total',
                Mortality != 0L) %>%
  dplyr::mutate(Log_population = log(Exposure),
                Zeros = ifelse(Age == 0, 1, 0),
                Age_fac = factor(Age),
                Sex = factor(Sex))

# View the observed age-mortality curves
ggplot(data = dat,
       aes(x = Age,
           y = Deaths / Exposure,
           group = Year,
           colour = Year)) +
  geom_line() +
  facet_wrap(~ Sex) +
  scale_colour_viridis_c() +
  labs(y = 'Observed log(Mortality)') +
  scale_y_log10()

# Fit a model to estimate how rates of mortality change over time
?fts
?ffc_gam
mod <- ffc_gam(
  formula = Deaths ~
    offset(log(Exposure)) +
    Sex * Zeros +
    # Use 10 thin-plate basis functions to characterize the age-mortality
    # function; for each of these functions, use a cubic
    # spline to estimate how it's coefficient changes through time
    fts(Age,
        bs = 'cr',
        k = 10,
        time_k = 15,
        time_bs = 'cr',
        time_m = 1,
        by = Sex),
  data = dat,
  time = 'Year',

  # Deaths is not an integer here, so use a Gamma distribution
  family = Gamma(link = 'log'),
  engine = 'bam',
  nthreads = 4,
  discrete = TRUE
)

# Have a look at the model object
summary(mod)
class(mod)
methods(class = 'ffc_gam')

# Draw some of the time-varying basis function coefficients
gratia::draw(
  mod,
  select = 1:2
)

# View predicted functional curves using a fixed offset
# to calculate a standardized rate of mortality
newdat <- dat
newdat$Exposure <- 1
newdat$preds <- predict(
  mod,
  newdata = newdat,
  type = 'response'
)

ggplot(data = newdat,
       aes(x = Age,
           y = preds,
           group = Year,
           colour = Year)) +
  geom_line() +
  facet_wrap(~ Sex) +
  scale_colour_viridis_c() +
  labs(y = 'Expected log10(Mortality)') +
  scale_y_log10()

# View the observed curves again as a sanity check
ggplot(data = dat,
       aes(x = Age,
           y = Mortality,
           group = Year,
           colour = Year)) +
  geom_line() +
  facet_wrap(~ Sex) +
  scale_colour_viridis_c() +
  labs(y = 'Observed log10(Mortality)') +
  scale_y_log10()

# Extract the time-varying coefficients, and their standard erros,
# in a tidy form
basis_ts <- fts_coefs(mod)
basis_ts

# To do:
# 1. Add functionality to plot the time-varying coefficients in a nice
#    and standard way
# 2. Add functionality to forecast the basis function coefficients, and to
#    combine these with the rest of the betas so that forecasts can be computed
# 3. Write unit tests and ensure the functions work with a broad range of bases,
#   including MRF smooths, cyclic smooths etc...
# 4. Update insight functions to ensure marginaleffects works without silly warnings
