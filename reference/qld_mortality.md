# Queensland mortality data

A dataset containing number of deaths, population at risk and ages for
residents of Queensland, Australia from 1979 to 2020

## Usage

``` r
qld_mortality
```

## Format

A `data.frame` containing the following fields:

- year:

  integer, year of records

- age:

  integer, age at death, in years

- sex:

  factor differentiating the sexes

- deaths:

  integer, number of deaths recorded

- population:

  numeric population recorded at 30 June each year

## Source

Australian Human Mortality Database

## Details

The age group 100 also includes people who died aged older than 100. The
data come from the Australian Human Mortality Database
(<https://aushd.org>).
