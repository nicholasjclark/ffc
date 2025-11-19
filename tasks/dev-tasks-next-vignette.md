# Development Task: El Niño Sea Surface Temperature Vignette

## Overview
Create a vignette demonstrating `ffc` package capabilities using El Niño sea surface temperature data from the `rainbow` package. This example showcases time-varying functional coefficients for seasonal patterns.

## Key Features to Demonstrate
- Time-varying functional coefficients using `fts()` terms
- Cyclic cubic splines for seasonal patterns
- GPDF forecasting with dynamic factor models
- Integration with `fable` ecosystem for forecast visualization
- Seasonal pattern analysis across multiple years

## Original Data Source
- **Package**: `rainbow::ElNino_OISST_region_1and2`
- **Structure**: Monthly SST data (12 months × 37 years)
- **Time Coverage**: 1982-2018
- **Variables**: Temperature measurements across months and years

## Model Specification

```r
mod_elnino <- ffc_gam(
  temperature ~
    fts(year, mean_only = TRUE, time_k = 10) +
    fts(month, bs = "cc", time_k = 25),
  data = train_tsibble,
  knots = list(month = c(0.5, 12.5)),
  time = "year",
  family = gaussian()
)
```

## Key Implementation Points

### Data Preparation
Create standalone R script to load data, process and save as .rda in data/
1. Convert matrix format to long-form tibble suitable for `ffc_gam`
2. Create proper time indexing for temporal dependencies
3. Split into training (≤2014) and test (>2014) sets
4. Convert to `tsibble` for `fable` integration
5. Save as .rda in data/ and add appropriate script in R/ to describe in roxygen2 documentation (vignette will work directly with the pre-loaded dataset)

### Model Features
- **Year effects**: `fts(year, mean_only = TRUE)` - captures inter-annual variability
- **Seasonal effects**: `fts(month, bs = "cc")` - cyclic cubic splines for seasonal patterns
- **Knots specification**: `knots = list(month = c(0.5, 12.5))` for proper cyclic boundaries

### Forecasting Approach
- **Method**: GPDF (Gaussian Process Dynamic Factor) model
- **Ensemble**: `mean_model = "ENS"` for improved accuracy
- **Factors**: `K = 3` latent factors
- **Validation**: Visual comparison of forecast vs actual seasonal curves

## Visualization Components

### 1. Seasonal Trend Analysis
- Monthly time series by year
- Color-coded by month to show seasonal patterns
- Comparison of actual vs forecast trajectories, showing forecast means and uncertainties

### 2. Annual Curve Comparison
- Seasonal curves grouped by year periods
- Overlay of actual and forecast annual patterns, showing means and uncertainties
- Assessment of forecast accuracy across different years

### 3. Model Diagnostics
- Coefficient evolution plots using `fts_coefs()`
- Model fit assessment with `gratia::draw()`
- Conditional predictions with `marginaleffects`

## Expected Outcomes

### Scientific Value
- Demonstrate climate pattern forecasting capabilities
- Show seasonal decomposition and inter-annual variability
- Validate model performance on real oceanographic data

### Package Demonstration
- Showcase advanced `fts()` term usage
- Illustrate cyclic spline functionality
- Demonstrate Stan-based dynamic factor forecasting
- Show integration with `tidyverse` and `fable` ecosystems

## Technical Considerations

### Performance
- Use `silent = TRUE` for Stan sampling to reduce output
- Limit iterations (`iter = 300`) for development testing
- Single chain (`chains = 1`) for faster prototyping

### Code Organization
- Modular structure suitable for vignette format
- Clear section breaks for explanation
- Comprehensive comments for educational value
- Error handling and data validation

## Implementation Status
- [x] Initial R script created with full workflow
- [ ] Convert to formal vignette format
- [ ] Add explanatory text and methodology
- [ ] Include interpretation of results
- [ ] Add references to El Niño/oceanography literature
- [ ] Test on different R configurations
- [ ] Integrate with pkgdown documentation

## Notes
This vignette would complement the existing mortality analysis vignette by:
- Showing environmental/climate applications
- Demonstrating cyclic spline functionality
- Using different forecasting models (GPDF vs others)
- Providing time series visualization examples
