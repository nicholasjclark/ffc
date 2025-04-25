# Forecasting from a ffc_gam object
library(ffc)
qld_train <- qld_mortality %>%
  dplyr::filter(year < 2016)
qld_test <- qld_mortality %>%
  dplyr::filter(year >= 2016)
mod <- ffc_gam(
  deaths ~
    offset(log(population)) +
    sex +
    fts(
      age,
      k = 10, bs = "cr", by = sex,
      time_bs = "cr", time_k = 10
    ),
  time = "year",
  data = qld_train,
  family = poisson(),
  engine = "bam"
)
summary(mod)

fc <- forecast(mod,
               newdata = qld_test)
fc
class(fc)
