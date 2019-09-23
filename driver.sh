


#######################
# Run 'evaluate_watershed.R'
#######################

# Run using Watershed approximate inference
model="Watershed_approximate"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
number_of_dimensions="3" # Can take on any real number greater than or equal to one
input_file="example_data/watershed_example_data.txt" # Input file
output_prefix="model_"$model"_number_of_dimensions_"$number_of_dimensions
Rscript evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model

# Run using Watershed exact inference
model="Watershed_exact"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
number_of_dimensions="3" # Can take on any real number greater than or equal to one
input_file="example_data/watershed_example_data.txt" # Input file
output_prefix="model_"$model"_number_of_dimensions_"$number_of_dimensions
Rscript evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model

# Run using RIVER
model="RIVER"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
number_of_dimensions="3" # Can take on any real number greater than or equal to one
input_file="example_data/watershed_example_data.txt" # Input file
output_prefix="model_"$model"_number_of_dimensions_"$number_of_dimensions
Rscript evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model




#######################
# Run 'predict_watershed.R'
# Note for convenience, the training file is the same as the prediction file. This does not necessarily have to be the case
#######################

# Run using Watershed approximate inference
model="Watershed_approximate"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
number_of_dimensions="3" # Can take on any real number greater than or equal to one
input_file="example_data/watershed_example_data.txt" # Input file
output_prefix="model_"$model"_number_of_dimensions_"$number_of_dimensions
Rscript predict_watershed.R --training_input $input_file --prediction_input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model

# Run using Watershed exact inference
model="Watershed_exact"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
number_of_dimensions="3" # Can take on any real number greater than or equal to one
input_file="example_data/watershed_example_data.txt" # Input file
output_prefix="model_"$model"_number_of_dimensions_"$number_of_dimensions
Rscript predict_watershed.R --training_input $input_file --prediction_input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model

# Run using RIVER
model="RIVER"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
number_of_dimensions="3" # Can take on any real number greater than or equal to one
input_file="example_data/watershed_example_data.txt" # Input file
output_prefix="model_"$model"_number_of_dimensions_"$number_of_dimensions
Rscript predict_watershed.R --training_input $input_file --prediction_input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model

