# Package index

## Core Functions

Main functions for functional forecasting

- [`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)
  : Fit a functional time series model using dynamic functional
  coefficients

- [`fts()`](https://nicholasjclark.github.io/ffc/reference/fts.md) :

  Define functions with dynamic coefficients in ffc formulae

## Model Methods

Methods for working with fitted models

- [`predict(`*`<ffc_gam>`*`)`](https://nicholasjclark.github.io/ffc/reference/predict.ffc_gam.md)
  :

  Predict from a fitted ffc_gam model

- [`forecast(`*`<ffc_gam>`*`)`](https://nicholasjclark.github.io/ffc/reference/forecast.ffc_gam.md)
  :

  Forecasting `ffc_gam` models

- [`forecast(`*`<fts_ts>`*`)`](https://nicholasjclark.github.io/ffc/reference/forecast.fts_ts.md)
  : Forecasting functional basis coefficients (Internal)

- [`model.frame(`*`<ffc_gam>`*`)`](https://nicholasjclark.github.io/ffc/reference/model.frame.ffc_gam.md)
  :

  Extract model.frame from a fitted `ffc_gam` object

## Coefficient Extraction

Extract and work with time-varying coefficients

- [`fts_coefs()`](https://nicholasjclark.github.io/ffc/reference/fts_coefs.ffc_gam.md)
  : Extract time-varying basis coefficients
- [`autoplot(`*`<fts_ts>`*`)`](https://nicholasjclark.github.io/ffc/reference/autoplot.fts_ts.md)
  : Plot time series of fts basis coefficients
- [`print(`*`<fts_ts>`*`)`](https://nicholasjclark.github.io/ffc/reference/print.fts_ts.md)
  : Print an fts_ts tibble

## Forecasting Models

Specialized forecasting models for time-varying coefficients

- [`ARDF()`](https://nicholasjclark.github.io/ffc/reference/ARDF.md) :
  Fit an autoregressive dynamic factor model
- [`GPDF()`](https://nicholasjclark.github.io/ffc/reference/GPDF.md) :
  Fit a GP dynamic factor model
- [`VARDF()`](https://nicholasjclark.github.io/ffc/reference/VARDF.md) :
  Fit a Vector autoregressive dynamic factor model

## Integration with fabletools

Convert ffc objects to work with the fabletools ecosystem

- [`as_fable(`*`<ffc_gam>`*`)`](https://nicholasjclark.github.io/ffc/reference/as_fable.ffc_gam.md)
  : Convert ffc_gam forecasts to fable object
- [`reexports`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`s`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`te`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`gam`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`bam`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`nb`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`negbin`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`betar`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`cnorm`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`ocat`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`scat`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`ziP`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`multinom`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`tw`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`ldTweedie`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`rTweedie`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`twlss`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`autoplot`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`forecast`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`generate`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  : Objects exported from other packages

## Data

Example datasets included with the package

- [`growth_data`](https://nicholasjclark.github.io/ffc/reference/growth_data.md)
  : Berkeley growth study data
- [`qld_mortality`](https://nicholasjclark.github.io/ffc/reference/qld_mortality.md)
  : Queensland mortality data
- [`elnino_sst`](https://nicholasjclark.github.io/ffc/reference/elnino_sst.md)
  : El Niño Sea Surface Temperature Data

## Re-exported Functions

Functions re-exported from mgcv and other packages for convenience

- [`reexports`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`s`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`te`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`gam`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`bam`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`nb`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`negbin`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`betar`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`cnorm`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`ocat`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`scat`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`ziP`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`multinom`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`tw`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`ldTweedie`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`rTweedie`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`twlss`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`autoplot`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`forecast`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  [`generate`](https://nicholasjclark.github.io/ffc/reference/reexports.md)
  : Objects exported from other packages
