library(ffc)
library(tsibble)
devtools::load_all()

# Create the same test data as the failing test
set.seed(42)
dat_grouped <- data.frame(
  y = c(rnorm(20, 0.5), rnorm(20, 0.8)),
  time = rep(1:20, 2),
  group = rep(c("A", "B"), each = 20)
)

# Fit the model
mod_grouped <- suppressWarnings(ffc_gam(
  formula = y ~ fts(time, bs = 'ad'),
  data = dat_grouped,
  time = "time",
  family = gaussian()
))

# Create test data
newdata_good <- data.frame(
  group = c("A", "A", "B", "B"),
  time = c(11, 12, 11, 12)
)

# Get intermediate objects for inspection
time_var <- "time"
intermed_coefs <- fts_coefs(mod_grouped, summary = FALSE, times = 25)

# Create functional_coefs like in the actual forecast function
shift_tail = function(x) x - tail(x, 1)
functional_coefs <- intermed_coefs %>%
  dplyr::group_by(.basis, .time) %>%
  dplyr::mutate(
    .mean = mean(.estimate)
  ) %>%
  dplyr::select(.basis, .time, !!time_var, .mean) %>%
  dplyr::distinct() %>%
  dplyr::ungroup() %>%
  dplyr::group_by(.basis) %>%
  dplyr::mutate(
    .mean = shift_tail(.mean)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::rename(.estimate = .mean) %>%
  structure(
    class = c("fts_ts", "tbl_df", "tbl", "data.frame"),
    time_var = attr(intermed_coefs, "time_var"),
    index = attr(intermed_coefs, "index"),
    index2 = attr(intermed_coefs, "index2")
  )

# Debug: Check the functional_coefs structure
print("=== functional_coefs structure ===")
print(class(functional_coefs))
print(head(functional_coefs))

# Debug: Do a small forecast to see what functional_fc looks like
functional_fc <- suppressWarnings(
  forecast(
    object = functional_coefs,
    h = 2,  # Just 2 time steps
    times = 25,
    model = "RW",
    stationary = FALSE
  )
)

print("=== functional_fc structure ===")
print(class(functional_fc))
print("Columns:", colnames(functional_fc))
print("Key vars:", key_vars(functional_fc))
print("Index:", index(functional_fc))
print("First few rows:")
print(head(functional_fc))

# Show what happens with time mapping
target_times <- newdata_good[[time_var]]  # [11, 12, 11, 12]
forecast_times <- sort(unique(functional_fc$.time))  # consecutive integers

print("=== Time mapping issue ===")
print("target_times:", target_times)
print("forecast_times:", forecast_times)
print("unique target_times:", unique(target_times))
print("unique forecast_times:", unique(forecast_times))

if(length(unique(forecast_times)) == length(unique(target_times))) {
  time_mapping <- setNames(unique(target_times), unique(forecast_times))
  print("time_mapping:", time_mapping)
  
  # Show what the remapped time column would look like
  remapped_times <- time_mapping[as.character(functional_fc$.time)]
  print("Original .time:", head(functional_fc$.time, 10))
  print("Remapped .time:", head(remapped_times, 10))
  
  # Check for duplicates in key+time combinations
  if(length(key_vars(functional_fc)) > 0) {
    key_cols <- key_vars(functional_fc)
    test_df <- functional_fc %>%
      as_tibble() %>%
      select(all_of(key_cols), .time) %>%
      mutate(.time = time_mapping[as.character(.time)])
    
    print("=== After remapping check ===")
    print("Key variables:", key_cols)
    duplicated_rows <- test_df %>%
      group_by(across(all_of(c(key_cols, ".time")))) %>%
      filter(n() > 1)
    
    if(nrow(duplicated_rows) > 0) {
      print("DUPLICATE key+time combinations found!")
      print(duplicated_rows)
    } else {
      print("No duplicates found in key+time combinations")
    }
  }
}