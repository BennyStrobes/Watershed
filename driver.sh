





if false; then

model="RIVER"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
number_of_dimensions="3" # Can take on any real number greater than or equal to one
input_file="example_data/watershed_example_data.txt" # Input file
pseudocount="10"
lambda="0.01" # Can take on a number or 'NA' to do a grid search
n2_pair_pvalue_fraction=".1"
binary_pvalue_thresh=".1"
output_prefix="model_"$model"_number_of_dimensions_"$number_of_dimensions"_lambda_"$lambda
Rscript evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model --dirichlet_prior_parameter $pseudocount --l2_prior_parameter $lambda --n2_pair_pvalue_fraction $n2_pair_pvalue_fraction --binary_pvalue_threshold $binary_pvalue_thresh

model="Watershed_approximate"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
output_prefix="model_"$model"_number_of_dimensions_"$number_of_dimensions"_lambda_"$lambda
Rscript evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model --dirichlet_prior_parameter $pseudocount --l2_prior_parameter $lambda --n2_pair_pvalue_fraction $n2_pair_pvalue_fraction --binary_pvalue_threshold $binary_pvalue_thresh

model="Watershed_exact"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
output_prefix="model_"$model"_number_of_dimensions_"$number_of_dimensions"_lambda_"$lambda
Rscript evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model --dirichlet_prior_parameter $pseudocount --l2_prior_parameter $lambda --n2_pair_pvalue_fraction $n2_pair_pvalue_fraction --binary_pvalue_threshold $binary_pvalue_thresh



model="RIVER"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
number_of_dimensions="1" # Can take on any real number greater than or equal to one
input_file="example_data/river_example_data_pheno_1.txt" # Input file
pseudocount="10"
lambda="0.01" # Can take on a number or 'NA' to do a grid search
n2_pair_pvalue_fraction=".1"
binary_pvalue_thresh=".1"
output_prefix="model_"$model"_number_of_dimensions_"$number_of_dimensions"_lambda_"$lambda"_pheno_1"
Rscript evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model --dirichlet_prior_parameter $pseudocount --l2_prior_parameter $lambda --n2_pair_pvalue_fraction $n2_pair_pvalue_fraction --binary_pvalue_threshold $binary_pvalue_thresh
fi



model="RIVER"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
number_of_dimensions="1" # Can take on any real number greater than or equal to one
input_file="example_data/river_example_data_pheno_2.txt" # Input file
pseudocount="10"
lambda="0.01" # Can take on a number or 'NA' to do a grid search
n2_pair_pvalue_fraction=".1"
binary_pvalue_thresh=".1"
output_prefix="model_"$model"_number_of_dimensions_"$number_of_dimensions"_lambda_"$lambda"_pheno_2"
Rscript evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model --dirichlet_prior_parameter $pseudocount --l2_prior_parameter $lambda --n2_pair_pvalue_fraction $n2_pair_pvalue_fraction --binary_pvalue_threshold $binary_pvalue_thresh


model="RIVER"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
number_of_dimensions="1" # Can take on any real number greater than or equal to one
input_file="example_data/river_example_data_pheno_3.txt" # Input file
pseudocount="10"
lambda="0.01" # Can take on a number or 'NA' to do a grid search
n2_pair_pvalue_fraction=".1"
binary_pvalue_thresh=".1"
output_prefix="model_"$model"_number_of_dimensions_"$number_of_dimensions"_lambda_"$lambda"_pheno_3"
Rscript evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model --dirichlet_prior_parameter $pseudocount --l2_prior_parameter $lambda --n2_pair_pvalue_fraction $n2_pair_pvalue_fraction --binary_pvalue_threshold $binary_pvalue_thresh
