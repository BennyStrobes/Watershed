#!/bin/R

# Test main functions in WatershedR 
library(WatershedR)

# For all examples, use example data that has 3 E outlier p-value columns,
# which corresponds to number_of_dimensions = 3
input = paste0("https://raw.githubusercontent.com/BennyStrobes/Watershed/",
     "master/example_data/watershed_example_data.txt")

# Run using Watershed approximate inference
res = evaluate_watershed(input_file = input,
  model_name = "Watershed_approximate",
  number_of_dimensions = 3,
  output_prefix = "watershed_approximate_n3")

# Run using Watershed exact inference
res = evaluate_watershed(input_file = input,
  model_name = "Watershed_exact",
  number_of_dimensions = 3,
  output_prefix = "watershed_exact_n3")

# Run using RIVER
res = evaluate_watershed(input_file = input,
  model_name = "RIVER",
  number_of_dimensions = 3,
  output_prefix = "river_n3")

# Run using Watershed approximate inference
res = predict_watershed(training_input = input,
  prediction_input = input,
  model_name = "Watershed_approximate",
  number_dimensions = 3,
  output_prefix = "watershed_approximate_n3")

# Run using Watershed exact inference
predict_watershed(training_input = input,
  prediction_input = input,
  model_name = "Watershed_exact",
  number_dimensions = 3,
  output_prefix = "watershed_exact_n3")

# Run using RIVER
predict_watershed(training_input = input,
  prediction_input = input,
  model_name = "RIVER",
  number_dimensions = 3,
  output_prefix = "river_n3")
