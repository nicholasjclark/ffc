devtools::load_all()
library(tsibble)
library(dplyr)

# Let me just create a simple tsibble similar to what the forecast function might create
# and see what happens when we try to remap the time values

# Create a mock tsibble with grouped data structure
mock_data <- data.frame(
  .time = c(1, 2, 1, 2),  # consecutive forecast times 
  .basis = c("basis1", "basis1", "basis2", "basis2"),
  .sim = c(0.1, 0.2, 0.3, 0.4)
)

# Convert to tsibble with .basis as key
mock_tsibble <- tsibble(
  mock_data,
  key = .basis,
  index = .time
)

print("=== Original mock_tsibble ===")
print(mock_tsibble)
print(paste("Key vars:", paste(key_vars(mock_tsibble), collapse = ", ")))
print(paste("Index:", index_var(mock_tsibble)))

# Now simulate the time mapping problem
# We want to map:
# 1 -> 11 (for both basis1 and basis2)  
# 2 -> 12 (for both basis1 and basis2)
# This creates duplicate (basis, time) combinations which should violate tsibble

target_times <- c(11, 12, 11, 12)  # From grouped data
forecast_times <- c(1, 2)  # From forecast

print("\n=== Time mapping ===")
print("target_times:", target_times)
print("forecast_times:", forecast_times)
print("unique(target_times):", unique(target_times))

# Create the mapping like in the actual code
time_mapping <- setNames(unique(target_times), forecast_times)
print("time_mapping:", time_mapping)

# Now try the mutate that's failing
print("\n=== Attempting mutate (should fail) ===")
tryCatch({
  result <- mock_tsibble %>%
    dplyr::mutate(.time = time_mapping[as.character(.time)])
  print("Success! Result:")
  print(result)
}, error = function(e) {
  print("Error occurred:")
  print(e$message)
})

# Show what the problem is - after remapping we have the same .time values
# for multiple .basis values, but the issue is that we're not accounting for
# the grouped nature of the data properly

print("\n=== Show what the remapped data would look like ===")
test_df <- as_tibble(mock_tsibble)
test_df$.time <- time_mapping[as.character(test_df$.time)]
print(test_df)

# Check for duplicates
duplicates <- test_df %>%
  group_by(.basis, .time) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)

if(nrow(duplicates) > 0) {
  print("Duplicate (key, time) combinations:")
  print(duplicates)
} else {
  print("No duplicate (key, time) combinations")
}