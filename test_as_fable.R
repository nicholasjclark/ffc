# Test script for as_fable.ffc_gam function with tourism data
# This script replicates the README example to debug the as_fable function

# Load packages
devtools::load_all()
library(fable)
library(tsibble)
library(dplyr)
library(ggplot2)

# Prepare the tourism data
tourism_melb <- tourism %>%
  filter(
    Region == "Melbourne",
    Purpose == "Visiting"
  ) %>%
  mutate(
    quarter = lubridate::quarter(Quarter),
    time = dplyr::row_number()
  )

print("Tourism data structure:")
str(tourism_melb)

# Split into training and testing folds
train <- tourism_melb %>%
  dplyr::slice_head(n = 75)

test <- tourism_melb %>%
  dplyr::slice_tail(n = 5)

print("Test data structure:")
str(test)

# Fit the ffc_gam model
mod <- ffc_gam(
  Trips ~
    # Use mean_only = TRUE to model a time-varying mean
    fts(
      time,
      mean_only = TRUE,
      time_k = 50,
      time_m = 1
    ) +
    # Time-varying seasonality
    fts(
      quarter,
      k = 4,
      time_k = 15,
      time_m = 1
    ),
  time = "time",
  data = train,
  family = tw(),
  engine = "gam"
)

print("Model summary:")
print(summary(mod))

# Generate forecasts
fc <- forecast(
  object = mod,
  newdata = test,
  model = "ETS",
  summary = FALSE
)

print("Forecast structure:")
str(fc)
print(class(fc))

# Try to convert to fable
print("Attempting as_fable conversion...")

# First check what columns test has
print("Test data columns:")
names(test)

# Check for time columns
print("Columns that are yearquarter/yearmonth:")
sapply(test, function(x) class(x))

# Try the conversion
tryCatch({
  fc_fable <- as_fable(mod, newdata = test, forecasts = fc)
  print("Success! Fable object created:")
  print(fc_fable)
  
  # Check fable structure
  print("Fable structure:")
  str(fc_fable)
  
  # Test that it works with autoplot
  print("Testing autoplot...")
  p <- autoplot(fc_fable, train)
  print("Autoplot successful")
  
  # Test that it works with accuracy scoring
  print("Testing accuracy scoring...")
  acc <- accuracy(fc_fable, test)
  print("Accuracy scoring successful:")
  print(acc)
  
  # Test other fabletools functions
  print("Testing other fabletools functions...")
  
  # Test forecast intervals
  print("Checking forecast intervals...")
  intervals <- hilo(fc_fable, level = c(80, 95))
  print("Forecast intervals successful:")
  print(intervals)
  
  # Test forecast combination (if multiple models)
  print("Testing forecast summary functions...")
  fc_summary <- fc_fable %>%
    summarise(
      mean_forecast = mean(.dist),
      median_forecast = median(.dist),
      q25 = quantile(.dist, 0.25),
      q75 = quantile(.dist, 0.75)
    )
  print("Forecast summary successful:")
  print(fc_summary)
  
}, error = function(e) {
  print(paste("Error:", e$message))
  print("Debugging info:")
  
  # Let's try to debug step by step
  print("Step 1: Check model time_var")
  print(paste("Model time_var:", mod$time_var))
  
  print("Step 2: Check if time_var exists in test")
  print(paste("time_var in test:", mod$time_var %in% names(test)))
  
  print("Step 3: Check response variable")
  response <- as.character(mod$formula[[2]])
  print(paste("Response:", response))
  print(paste("Response in test:", response %in% names(test)))
  
  print("Step 4: Try manual fable construction")
  # Try to manually construct the fable
  fable_data <- test
  fable_data[[".dist"]] <- fc
  fable_data[[".mean"]] <- as.numeric(mean(fc))
  
  print("Manual fable_data structure before tsibble:")
  str(fable_data)
  
  # Try tsibble conversion
  tryCatch({
    # Check what time index to use
    if ("Quarter" %in% names(fable_data) && inherits(fable_data$Quarter, "yearquarter")) {
      print("Using Quarter as index")
      fable_tsbl <- tsibble::as_tsibble(
        fable_data,
        index = Quarter
      )
    } else if ("time" %in% names(fable_data)) {
      print("Using time as index")
      fable_tsbl <- tsibble::as_tsibble(
        fable_data,
        index = time
      )
    }
    
    print("Tsibble created successfully:")
    print(fable_tsbl)
    
    # Now try build_fable
    print("Trying build_fable...")
    fc_fable_manual <- fabletools:::build_fable(
      fable_tsbl,
      response = "Trips",
      distribution = ".dist"
    )
    
    print("Manual fable created successfully!")
    print(fc_fable_manual)
    
  }, error = function(e2) {
    print(paste("Error in manual construction:", e2$message))
  })
})