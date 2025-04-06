library(vital)
library(tidyverse); theme_set(theme_classic(base_size = 12))
library(mgcv)
library(marginaleffects)

# Experimenting with dynamic factor model for functional curve
# forecasting
dat <- aus_mortality

dat <- dat %>%
  dplyr::filter(Year > 1970,
                State == 'Queensland',
                Sex != 'total') %>%
  dplyr::mutate(Log_population = log(Exposure),
                Zeros = ifelse(Age == 0, 1, 0),
                Deaths = ceiling(Deaths),
                Sex = factor(Sex))

# Create a set of basis functions for the target (Age)
base_smooth <- smoothCon(
  s(Age, bs = 'tp', k = 11),
  data = data.frame(
    Age = unique(dat$Age)
  ),
  null.space.penalty = TRUE
)[[1]]

# Now create the design matrix of basis functions
X <- as.data.frame(
  PredictMat(
    base_smooth,
    data = dat
  )
)
nulls <- apply(X, 2, sd) == 0
X <- X[, !nulls]
colnames(X) <- paste0('bfun_',
                      1:NCOL(X),
                      '_male')

# Bind to original data
dat <- bind_cols(
  dat, X
)

# Now add the female functions
colnames(X) <- paste0('bfun_',
                      1:NCOL(X),
                      '_female')
dat <- bind_cols(
  dat, X
)

# Build the hierarchical GAM formula
# where male and female time series share
# smoothing parameters
mod <- bam(
  Deaths ~
    offset(Log_population) +
    Sex * Zeros +
    Age +
    s(Year, by = bfun_1_male, id = 1, bs = 'cr') +
    s(Year, by = bfun_2_male, id = 2, bs = 'cr') +
    s(Year, by = bfun_3_male, id = 3, bs = 'cr') +
    s(Year, by = bfun_4_male, id = 4, bs = 'cr') +
    s(Year, by = bfun_5_male, id = 5, bs = 'cr') +
    s(Year, by = bfun_6_male, id = 6, bs = 'cr') +
    s(Year, by = bfun_7_male, id = 7, bs = 'cr') +
    s(Year, by = bfun_8_male, id = 8, bs = 'cr') +
    s(Year, by = bfun_9_male, id = 9, bs = 'cr') +
    s(Year, by = bfun_10_male, id = 10, bs = 'cr') +

    s(Year, by = bfun_1_female, id = 1, bs = 'cr') +
    s(Year, by = bfun_2_female, id = 2, bs = 'cr') +
    s(Year, by = bfun_3_female, id = 3, bs = 'cr') +
    s(Year, by = bfun_4_female, id = 4, bs = 'cr') +
    s(Year, by = bfun_5_female, id = 5, bs = 'cr') +
    s(Year, by = bfun_6_female, id = 6, bs = 'cr') +
    s(Year, by = bfun_7_female, id = 7, bs = 'cr') +
    s(Year, by = bfun_8_female, id = 8, bs = 'cr') +
    s(Year, by = bfun_9_female, id = 9, bs = 'cr') +
    s(Year, by = bfun_10_female, id = 10, bs = 'cr'),
  data = dat,
  family = poisson(),
  discrete = TRUE
)
summary(mod)
gratia::draw(
  mod,
  select = c(1:3,
             11:13)
)

plot_predictions(
  mod,
  by = c('Age', 'Year', 'Sex'),
  type = 'link'
) +
  labs(y = 'log(mortality) rate',
       title = 'Standardised mortality rates')

plot_predictions(
  mod,
  by = c('Year', 'bfun_3_male'),
  newdata = datagrid(
    bfun_3_male = min,
    Year = 1990:2020
  ),
  type = 'link'
)

# Extract basis coefficients for forecasting by using
# gratia; then just need to multiply the forecasted coefs
# by a matrix of 1s
names(coef(mod))
coef_idx <- grepl('Year', names(coef(mod)))

lpmat <- predict(
  mod,
  type = 'lpmatrix'
)
length(which(lpmat[,181] != 0))
