# Berkeley growth study data

A dataset containing the heights of 39 boys measured at ages 3 to 15
years with regular yearly intervals. This is a subset of the original
Berkeley Growth Study data filtered to contain only regular time
intervals suitable for time series modeling.

## Usage

``` r
growth_data
```

## Format

A `data.frame` with 507 observations and 3 variables:

- height_cm:

  numeric height, in centimeters

- age_yr:

  numeric age at which the boys were measured, in years (3, 4, 5, ...,
  15)

- id:

  factor indicating the different subjects (39 boys)

## References

Tuddenham, R. D., and Snyder, M. M. (1954) "Physical growth of
California boys and girls from birth to age 18", *University of
California Publications in Child Development*, 1, 183-364.
