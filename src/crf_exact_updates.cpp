#include <vector>
#include <math.h>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <time.h>
#include <Rcpp.h>
using namespace Rcpp;


// Create matrix of dimension 2^(number_of_dimensions) X number_of_dimensions
// Each row of this matrix summarizes all possible binary values that CRF can take on
// [[Rcpp::export]]
NumericMatrix extract_all_binary_combinations(int n) {
	// Initialize output matrix
	int nTemp = (int)pow(2, n) - 1;
	NumericMatrix combo_mat(nTemp + 1, n);

	// Loop through all possible values the CRF can take on 
	for (int i = 0; i <= nTemp; i++) {
		// Loop through dimensions
		for (int k = 0; k < n; k++) {
			if ((i >> k) & 0x1){
				combo_mat(i,k) = 1;
			} else {
				combo_mat(i,k) = 0;
			}
		}
	}
	return combo_mat;
}

// Create matrix of dimension 2^(number_of_dimensions-1) X number_of_dimensions
// Each row of this matrix summarizes all possible binary values that CRF can take 
// on assuming the values of in column column_to_ignore are fixed to be 1
// [[Rcpp::export]]
NumericMatrix extract_all_binary_combinations_ignoring_one_row(int n, int column_to_ignore) {
	// Ignoring column_to_ignore, get matrix of dimension 2^(number_of_dimensions-1) X (number_of_dimensions -1) 
	// that representing all possible values CRF can take on for remaining dimensions
	NumericMatrix all_binary_combinations_matrix = extract_all_binary_combinations(n-1);

	// Initialize output matrix
	int N = all_binary_combinations_matrix.nrow();
	NumericMatrix combo_mat(N, n);
	// Loop through all possible values the CRF can take on 
	for (int row_num = 0; row_num < N; row_num++) {
		int counter = 0;
		// Loop through each dimension
		for (int column_num = 0; column_num < n; column_num++) {
			// If the dimension is the column_to_ignore, set it eqaul to one
			if (column_num == column_to_ignore) {
				combo_mat(row_num, column_num) = 1;
			// If the column is not a column to ignore (ie a dimension we are marginalizing out)
			} else {
				combo_mat(row_num, column_num) = all_binary_combinations_matrix(row_num, counter);
				counter += 1;
			}
		}
	}
	return combo_mat;
}

// Create matrix of dimension 2^(number_of_dimensions-2) X number_of_dimensions
// Each row of this matrix summarizes all possible binary values that CRF can take 
// on assuming the values of column_to_ignore1 and column_to_ignore2 are fixed to be 1
// [[Rcpp::export]]
NumericMatrix extract_all_binary_combinations_ignoring_two_row(int n, int column_to_ignore1, int column_to_ignore2) {
	// Special case (handeled seperately) for when there are only two dimensions
	if (n == 2) {
		NumericMatrix combo_mat(1,2);
		combo_mat(0,0) = 1;
		combo_mat(0,1) = 1;
		return combo_mat;
	// General case for when the number of dimensions is greater than 2
	} else {
		// Ignoring column_to_ignore1 and column_to_ignore2, get matrix of dimension 2^(number_of_dimensions-2) X (number_of_dimensions -2) that representing all possible values CRF can take on for remaining dimensions
		NumericMatrix all_binary_combinations_matrix = extract_all_binary_combinations(n-2);
		// Initialize output matrix
		int nrow = all_binary_combinations_matrix.nrow();
		NumericMatrix combo_mat(nrow, n);
		// Loop through all possible values CRF can take on
		for (int row_num = 0; row_num < nrow; row_num++) {
			int counter = 0;
			// Loop through each dimension
			for (int column_num = 0; column_num < n; column_num++) {
				// If the dimension is a column to ignore, set it equal to one
				if (column_num == column_to_ignore1 || column_num == column_to_ignore2) {
					combo_mat(row_num, column_num) = 1;
				// If the column is not a column to ignore (ie a dimension we are marginalizing out)
				} else {
					combo_mat(row_num, column_num) = all_binary_combinations_matrix(row_num, counter);
					counter += 1;
				}
			}
		}
		return combo_mat;
	}
}

// For a CRF value (particular row in all_binary_combinations_matrix), compute the relative CRF weight
// [[Rcpp::export]]
double un_normalized_crf_weight(NumericMatrix all_binary_combinations_matrix, int combination_number, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, bool posterior_bool) {
	// Initialize weight
	double weight = 0;
	int dimension_counter = 0;
	// Loop through dimensions
	for (int dimension=0; dimension < number_of_dimensions; dimension++) {
		// Add term from intercept
		weight += all_binary_combinations_matrix(combination_number, dimension)*theta_singleton(dimension);
		// Loop through features
		for (int d = 0; d < feat.ncol(); d++) {
			// Add term from features
			weight += all_binary_combinations_matrix(combination_number,dimension)*feat(sample_num,d)*theta(d,dimension);
		}
		// Loop through all pairs of dimensions
		for (int dimension2=dimension; dimension2 < number_of_dimensions; dimension2++) {
			if (dimension != dimension2) {
				// Add edge weight
				weight += all_binary_combinations_matrix(combination_number,dimension)*all_binary_combinations_matrix(combination_number,dimension2)*theta_pair(0, dimension_counter);
				dimension_counter += 1;
			}
		}
		// Check to see if we are supposed to incorperate expression data && whether the expression data is observed
		if (posterior_bool == true && discrete_outliers(sample_num, dimension) == discrete_outliers(sample_num, dimension)) {
			if (all_binary_combinations_matrix(combination_number, dimension) == 1) {
				weight += log(phi_outlier(dimension, discrete_outliers(sample_num, dimension) - 1));
			} else {
				weight += log(phi_inlier(dimension, discrete_outliers(sample_num, dimension) - 1));
			}
		}
	}
	return weight;
}

// Compute CRF normalization constant for a specifc sample
// [[Rcpp::export]]
double exact_normalization_constant(NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, bool posterior_bool) {
	// Extract matrix summarizing all possible values the CRF can take on
	// Create matrix of dimension 2^(number_of_dimensions) X number_of_dimensions
	// Each row of this matrix summarizes all possible binary values that CRF can take on
	NumericMatrix all_binary_combinations_matrix = extract_all_binary_combinations(number_of_dimensions);
	// Initialize variable keep track of normalization constant
	double val = 0;
	// Loop through each possible value the CRF can take on 
	for (int combination_number = 0; combination_number < all_binary_combinations_matrix.nrow(); combination_number++) {
		// And compute un-normalized CRF weight corresponding to that sample
		double un_normalized_wight = un_normalized_crf_weight(all_binary_combinations_matrix, combination_number, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, posterior_bool);
		val += exp(un_normalized_wight);
	}
	return log(val);
}

// Compute probability of a CRF label (Z*) 
// ie compute P(Z=Z*)
// [[Rcpp::export]]
double exact_probability(double normalization_constant, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, int combination_number,NumericMatrix all_binary_combinations_ignoring_one_row, bool posterior_bool) {
	double prob = exp(un_normalized_crf_weight(all_binary_combinations_ignoring_one_row, combination_number, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, posterior_bool) - normalization_constant);
	return prob;
}

// Compute marginal posterior probability for this sample in this dimension
// Compute P(Z_dimension=1|G) if posterior_bool==false
// Compute P(Z_dimension=1|G,E) if posterior_bool==true
// Involves marginalizing out all other dimensions
// [[Rcpp::export]]
double exact_marginal_probability(double normalization_constant, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, int dimension, bool posterior_bool) {
	// Create matrix of dimension 2^(number_of_dimensions-1) X number_of_dimensions
	// Each row of this matrix summarizes all possible binary values that CRF can take on assuming the values of in column column_to_ignore are fixed to be 1
	NumericMatrix all_binary_combinations_ignoring_one_row = extract_all_binary_combinations_ignoring_one_row(number_of_dimensions, dimension);
	double marginal_prob = 0;
	// Loop through each of the rows (possible CRF values) of the above matrix
	for (int combination_number = 0; combination_number < all_binary_combinations_ignoring_one_row.nrow(); combination_number++) {
		// For this CRF value, compute the probability
		marginal_prob += exact_probability(normalization_constant, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, combination_number, all_binary_combinations_ignoring_one_row, posterior_bool);
	}
	return marginal_prob;
}


// Compute marginal pairwise probability for dimension pair for a particular sample
// Compute P(Z_dimension1=1, Z_dimension2=1| G) if posterior_bool==false
// Compute P(Z_dimension1=1, Z_dimension2=1| G,E) if posterior_bool==true
// Involves marginalizing out all other dimensions
// [[Rcpp::export]]
double exact_marginal_pairwise_probability(double normalization_constant, int dimension1, int dimension2, int dimension_counter, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, bool posterior_bool) {
	// Initialize marginal probability variable
	double marginal_prob = 0;
	// Create matrix of dimension 2^(number_of_dimensions-2) X number_of_dimensions
	// Each row of this matrix summarizes all possible binary values that CRF can take on assuming the values of dimension_1 and dimension_2 are fixed to be 1
	NumericMatrix marginal_binary_combinations_matrix = extract_all_binary_combinations_ignoring_two_row(number_of_dimensions, dimension1, dimension2);
	// Loop through each of the rows (possible CRF values) of the above matrix
	for (int combination_number = 0; combination_number < marginal_binary_combinations_matrix.nrow(); combination_number++) {
		// For this CRF value, compute the probability
		marginal_prob += exact_probability(normalization_constant, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, combination_number, marginal_binary_combinations_matrix, posterior_bool);
	}
	return marginal_prob;
}

// Compute both marginal Posterior probability as well as marginal pairwise probabilities for Conditional random field using exact inference
// If posterior_bool==true, compute P(Z|E,G)
// If posterior_bool==false, compute P(Z|G)
// [[Rcpp::export]]
List update_marginal_probabilities_exact_inference_cpp(NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int number_of_pairs, bool posterior_bool) {
	// Initialize output matrices
	NumericMatrix probabilities(feat.nrow(), number_of_dimensions);
	NumericMatrix probabilities_pairwise(feat.nrow(), number_of_pairs);

	// Loop through samples
	for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
		// Compute normalization constant for this sample
		double normalization_constant = exact_normalization_constant(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, posterior_bool);

		// Initialize variable to help when iterating through pairs of edges
		int dimension_counter = 0;
		// Loop through dimensions
		for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
			// Compute marginal posterior probability for this sample in this dimension
			probabilities(sample_num, dimension) = exact_marginal_probability(normalization_constant, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, dimension, posterior_bool);
			// Loop through pairs of dimensions
			for (int dimension2=dimension; dimension2 < number_of_dimensions; dimension2++) {
				if (dimension != dimension2) {
					// Compute marginal pairwise probability for dimension pair
					probabilities_pairwise(sample_num, dimension_counter) = exact_marginal_pairwise_probability(normalization_constant, dimension, dimension2, dimension_counter, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, posterior_bool);
					dimension_counter += 1;
				}
			}
		}
	}
	List ret;
	ret["probability"] = probabilities;
	ret["probability_pairwise"] = probabilities_pairwise;
	return ret;
}



// Compute exact likelihood of K=number_of_dimensions dimensionsal Conditional Random Field (CRF)
// [[Rcpp::export]]
double compute_crf_likelihood_exact_inference_cpp(NumericMatrix posterior, NumericMatrix posterior_pairwise, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, double lambda, double lambda_pair, double lambda_singleton) {
	// Initialize output likelihood
	double log_likelihood = 0;
	// Loop through samples
	for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
		// Compute normalization constant for this sample
		double normalization_constant = exact_normalization_constant(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, false);
		// Subtract normalization constant from the log likelihood
		log_likelihood = log_likelihood - normalization_constant;
		// Loop through dimensions
		int dimension = 0;
		int dimension_counter = 0;
		for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
			// Add contribution of intercept term for this dimension
			log_likelihood += theta_singleton(dimension)*posterior(sample_num, dimension);
			// Add contribution of feature terms for this dimension
			for (int d = 0; d < feat.ncol(); d++) {
				log_likelihood += theta(d, dimension)*feat(sample_num, d)*posterior(sample_num, dimension);
			}
			// This nested for loop will loop through all pairs of dimensions
			for (int dimension2=dimension; dimension2 < number_of_dimensions; dimension2++) {
				if (dimension != dimension2) {
					// Add contribution of a given edge pair
					log_likelihood += theta_pair(0, dimension_counter)*posterior_pairwise(sample_num, dimension_counter);
					dimension_counter += 1;
				}
			}
		}
	}

	// Normalize log likelihood by the number of samples
	log_likelihood = log_likelihood/feat.nrow();

	// Add L2 penalties
	int dimension_counter = 0;
	for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
		// Generally do not regularize intercept terms (ie lambda_singleton is usually set to zero)
		log_likelihood = log_likelihood - .5*lambda_singleton*(theta_singleton(dimension)*theta_singleton(dimension));
		// Regularize edge weights
		for (int dimension2=dimension; dimension2 < number_of_dimensions; dimension2++) {
			if (dimension != dimension2) {
				log_likelihood = log_likelihood - .5*lambda_pair*(theta_pair(0,dimension_counter)*theta_pair(0, dimension_counter));
				dimension_counter += 1;
			}
		}
		// Regularize feature vectors 
		for (int d = 0; d < feat.ncol(); d++) {
			log_likelihood = log_likelihood - .5*lambda*(theta(d,dimension)*theta(d,dimension));
		}
	}
	return log_likelihood;
}
