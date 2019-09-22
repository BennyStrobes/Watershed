#include <vector>
#include <math.h>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <time.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix extract_all_binary_combinations(int n) {
	// Initialize output
	int nTemp = (int)pow(2, n) - 1;

	NumericMatrix combo_mat(nTemp + 1, n);


	for (int i = 0; i <= nTemp; i++)
	{
		for (int k = 0; k < n; k++)
		{
			if ((i >> k) & 0x1)
			{
				combo_mat(i,k) = 1;
			}
			else
			{
				combo_mat(i,k) = 0;
			}
		}
	}
	return combo_mat;
}

// NumericMatrix all_binary_combinations_ignoring_one_row = extract_all_binary_combinations_ignoring_one_row(4,3);
NumericMatrix extract_all_binary_combinations_ignoring_one_row(int n, int column_to_ignore) {
	NumericMatrix all_binary_combinations_matrix = extract_all_binary_combinations(n-1);

	int N = all_binary_combinations_matrix.nrow();
	NumericMatrix combo_mat(N, n);
	for (int row_num = 0; row_num < N; row_num++) {
		int counter = 0;
		for (int column_num = 0; column_num < n; column_num++) {
			if (column_num == column_to_ignore) {
				combo_mat(row_num, column_num) = 1;
			} else {
				combo_mat(row_num, column_num) = all_binary_combinations_matrix(row_num, counter);
				counter += 1;
			}
		}
	}
	return combo_mat;
}

// NumericMatrix all_binary_combinations_ignoring_two_row = extract_all_binary_combinations_ignoring_two_row(2, 1, 0);
NumericMatrix extract_all_binary_combinations_ignoring_two_row(int n, int column_to_ignore1, int column_to_ignore2) {
	if (n == 2) {
		NumericMatrix combo_mat(1,2);
		combo_mat(0,0) = 1;
		combo_mat(0,1) = 1;
		return combo_mat;
	} else {
		NumericMatrix all_binary_combinations_matrix = extract_all_binary_combinations(n-2);
		int nrow = all_binary_combinations_matrix.nrow();
		NumericMatrix combo_mat(nrow, n);
		for (int row_num = 0; row_num < nrow; row_num++) {
			int counter = 0;
			for (int column_num = 0; column_num < n; column_num++) {
				if (column_num == column_to_ignore1 || column_num == column_to_ignore2) {
					combo_mat(row_num, column_num) = 1;
				} else {
					combo_mat(row_num, column_num) = all_binary_combinations_matrix(row_num, counter);
					counter += 1;
				}
			}
		}
		return combo_mat;
	}
}

double un_normalized_crf_weight(NumericMatrix all_binary_combinations_matrix, int combination_number, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, bool posterior_bool) {
	// Initialize weight
	double weight = 0;
	int dimension_counter = 0;
	for (int dimension=0; dimension < number_of_dimensions; dimension++) {
		weight += all_binary_combinations_matrix(combination_number, dimension)*theta_singleton(dimension);
		// Loop through features
		for (int d = 0; d < feat.ncol(); d++) {
			weight += all_binary_combinations_matrix(combination_number,dimension)*feat(sample_num,d)*theta(d,dimension);
		}
		for (int dimension2=dimension; dimension2 < number_of_dimensions; dimension2++) {
			if (dimension != dimension2) {
				for (int theta_pair_dimension=0; theta_pair_dimension < theta_pair.nrow(); theta_pair_dimension++) {
					if (theta_pair_dimension == 0) {
						weight += all_binary_combinations_matrix(combination_number,dimension)*all_binary_combinations_matrix(combination_number,dimension2)*theta_pair(theta_pair_dimension, dimension_counter);
					} else {
						weight += feat(sample_num,(theta_pair_dimension-1))*all_binary_combinations_matrix(combination_number,dimension)*all_binary_combinations_matrix(combination_number,dimension2)*theta_pair(theta_pair_dimension, dimension_counter);
					}
				}
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

double exact_normalization_constant(NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, bool posterior_bool) {
	NumericMatrix all_binary_combinations_matrix = extract_all_binary_combinations(number_of_dimensions);
	double val = 0;

	for (int combination_number = 0; combination_number < all_binary_combinations_matrix.nrow(); combination_number++) {
		double un_normalized_wight = un_normalized_crf_weight(all_binary_combinations_matrix, combination_number, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, posterior_bool);
		val += exp(un_normalized_wight);
	}
	return log(val);
}

double exact_probability(double normalization_constant, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, int combination_number,NumericMatrix all_binary_combinations_ignoring_one_row, bool posterior_bool) {
	double prob = exp(un_normalized_crf_weight(all_binary_combinations_ignoring_one_row, combination_number, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, posterior_bool) - normalization_constant);
	return prob;
}

double exact_marginal_probability(double normalization_constant, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, int dimension, bool posterior_bool) {
	NumericMatrix all_binary_combinations_ignoring_one_row = extract_all_binary_combinations_ignoring_one_row(number_of_dimensions, dimension);
	double marginal_prob = 0;
	for (int combination_number = 0; combination_number < all_binary_combinations_ignoring_one_row.nrow(); combination_number++) {
		marginal_prob += exact_probability(normalization_constant, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, combination_number, all_binary_combinations_ignoring_one_row, posterior_bool);
	}
	return marginal_prob;
}


double exact_observed_sample_likelihood(double normalization_constant, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num) {
	NumericMatrix all_binary_combinations_matrix = extract_all_binary_combinations(number_of_dimensions);
	double prob = 0;
	for (int combination_number = 0; combination_number < all_binary_combinations_matrix.nrow(); combination_number++) {
		double combination_prob = exact_probability(normalization_constant, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, combination_number, all_binary_combinations_matrix, false);
		for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
			// Check to make sure expression is measured
			if (discrete_outliers(sample_num, dimension) == discrete_outliers(sample_num, dimension)) {
				if (all_binary_combinations_matrix(combination_number, dimension) == 1) {
					combination_prob = combination_prob*phi_outlier(dimension, discrete_outliers(sample_num, dimension) - 1);
				} else {
					combination_prob = combination_prob*phi_inlier(dimension, discrete_outliers(sample_num, dimension) - 1);
				}
			}
		}
		prob += combination_prob;
	}
	return prob;
}


double exact_marginal_pairwise_probability(double normalization_constant, int dimension1, int dimension2, int dimension_counter, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, bool posterior_bool) {
	double marginal_prob = 0;
	NumericMatrix marginal_binary_combinations_matrix = extract_all_binary_combinations_ignoring_two_row(number_of_dimensions, dimension1, dimension2);
	for (int combination_number = 0; combination_number < marginal_binary_combinations_matrix.nrow(); combination_number++) {
		marginal_prob += exact_probability(normalization_constant, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, combination_number, marginal_binary_combinations_matrix, posterior_bool);
	}

	return marginal_prob;
}

// [[Rcpp::export]]
List update_marginal_probabilities_exact_inference_cpp(NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int number_of_pairs, bool posterior_bool) {
	// Rcpp::Rcout << discrete_outliers(0,1) << std::endl;  
	// bool temp = discrete_outliers(0,1) == discrete_outliers(0,1);
	// Rcpp::Rcout << temp << std::endl;  

	// Initialize output matrices
	NumericMatrix probabilities(feat.nrow(), number_of_dimensions);
	NumericMatrix probabilities_pairwise(feat.nrow(), number_of_pairs);

	// Loop through samples
	for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
		// Compute normalization constant for this sample
		// normalization_constant <- exact_posterior_normalization_constant(feat[n,], discrete_outliers[n,], model_params$theta_singleton, model_params$theta_pair, model_params$theta, model_params$phi, model_params$number_of_dimensions)
		double normalization_constant = exact_normalization_constant(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, posterior_bool);
		// Loop through dimensions
		int dimension = 0;
		int dimension_counter = 0;
		for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
			probabilities(sample_num, dimension) = exact_marginal_probability(normalization_constant, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, dimension, posterior_bool);
			for (int dimension2=dimension; dimension2 < number_of_dimensions; dimension2++) {
				if (dimension != dimension2) {
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

// [[Rcpp::export]]
List compute_all_exact_posterior_predictions_cpp(NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions) {
	// All output combinations
	NumericMatrix combo_mat = extract_all_binary_combinations(number_of_dimensions);
	// Initialize output matrices
	NumericMatrix probabilities(feat.nrow(), combo_mat.nrow());
	// Loop through samples
	for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
		// Compute normalization constant for this sample
		double normalization_constant = exact_normalization_constant(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, true);
		// Loop through all possible states
		for (int combo_num = 0; combo_num < combo_mat.nrow(); combo_num++) {
			probabilities(sample_num, combo_num) = exact_probability(normalization_constant, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, combo_num, combo_mat, true);
		}
	}
	List ret;
	ret["probability"] = probabilities;
	ret["combination"] = combo_mat;
	return ret;
}

// [[Rcpp::export]]
List compute_all_exact_crf_posterior_predictions_cpp(NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions) {
	// All output combinations
	NumericMatrix combo_mat = extract_all_binary_combinations(number_of_dimensions);
	// Initialize output matrices
	NumericMatrix probabilities(feat.nrow(), combo_mat.nrow());
	// Loop through samples
	for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
		// Compute normalization constant for this sample
		double normalization_constant = exact_normalization_constant(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, false);
		// Loop through all possible states
		for (int combo_num = 0; combo_num < combo_mat.nrow(); combo_num++) {
			probabilities(sample_num, combo_num) = exact_probability(normalization_constant, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, combo_num, combo_mat, false);
		}
	}
	List ret;
	ret["probability"] = probabilities;
	ret["combination"] = combo_mat;
	return ret;
}

// [[Rcpp::export]]
double compute_exact_observed_data_log_likelihood_cpp(NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int number_of_pairs) {
	// Initialize log_likelihood
	double log_likelihood = 0;
	// Loop through samples
	for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
		// Compute normalization constant for this sample
		double normalization_constant = exact_normalization_constant(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, false);
		// Compute observed data likelihood for this sample
		double observed_sample_likelihood = exact_observed_sample_likelihood(normalization_constant, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num);
		// Add log of this observed data likelihood to a variable log_likelihood
		log_likelihood += log(observed_sample_likelihood);
	}
	// log_likelihood = log_likelihood/feat.nrow();
	return log_likelihood;
}

// [[Rcpp::export]]
double compute_crf_likelihood_exact_inference_cpp(NumericMatrix posterior, NumericMatrix posterior_pairwise, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, double lambda, double lambda_pair, double lambda_singleton) {
	// Rcpp::Rcout << feat(1,1) << std::endl;  
	// Initialize output likelihood
	double log_likelihood = 0;
	// Loop through samples
	for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
		// Compute normalization constant for this sample
		double normalization_constant = exact_normalization_constant(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, false);
		log_likelihood = log_likelihood - normalization_constant;
		// Loop through dimensions
		int dimension = 0;
		int dimension_counter = 0;
		for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
			log_likelihood += theta_singleton(dimension)*posterior(sample_num, dimension);
			for (int d = 0; d < feat.ncol(); d++) {
				log_likelihood += theta(d, dimension)*feat(sample_num, d)*posterior(sample_num, dimension);
			}
			for (int dimension2=dimension; dimension2 < number_of_dimensions; dimension2++) {
				if (dimension != dimension2) {
					for (int theta_pair_dimension=0; theta_pair_dimension < theta_pair.nrow(); theta_pair_dimension++) {
						if (theta_pair_dimension == 0) {
							log_likelihood += theta_pair(theta_pair_dimension, dimension_counter)*posterior_pairwise(sample_num, dimension_counter);
						} else {
							log_likelihood += feat(sample_num, (theta_pair_dimension-1))*theta_pair(theta_pair_dimension, dimension_counter)*posterior_pairwise(sample_num, dimension_counter);
						}
					}
					dimension_counter += 1;
				}
			}
		}
	}

	log_likelihood = log_likelihood/feat.nrow();

	// Add L2 penalties
	int dimension_counter = 0;
	for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
		log_likelihood = log_likelihood - .5*lambda_singleton*(theta_singleton(dimension)*theta_singleton(dimension));
		for (int dimension2=dimension; dimension2 < number_of_dimensions; dimension2++) {
			if (dimension != dimension2) {
				log_likelihood = log_likelihood - .5*lambda_pair*(theta_pair(0,dimension_counter)*theta_pair(0, dimension_counter));
				for (int theta_pair_dimension=1; theta_pair_dimension < theta_pair.nrow(); theta_pair_dimension++) {
					log_likelihood = log_likelihood - .5*lambda*(theta_pair(theta_pair_dimension,dimension_counter)*theta_pair(theta_pair_dimension, dimension_counter));
				}
				dimension_counter += 1;
			}
		}
		for (int d = 0; d < feat.ncol(); d++) {
			log_likelihood = log_likelihood - .5*lambda*(theta(d,dimension)*theta(d,dimension));
		}
	}
	return log_likelihood;
}
