#include <vector>
#include <math.h>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <time.h>
#include <Rcpp.h>
using namespace Rcpp;


double un_normalized_independent_crf_weight(int dimension, int combination_number, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, bool posterior_bool) {
	// Initialize weight
	double weight = 0;
	int dimension_counter = 0;
	weight += combination_number*theta_singleton(dimension);
	for (int d = 0; d < feat.ncol(); d++) {
		weight += combination_number*feat(sample_num,d)*theta(d,dimension);
	}
	if (posterior_bool == true && discrete_outliers(sample_num, dimension) == discrete_outliers(sample_num, dimension)) {
		if (combination_number == 1) {
			weight += log(phi_outlier(dimension, discrete_outliers(sample_num, dimension) - 1));
		} else {
			weight += log(phi_inlier(dimension, discrete_outliers(sample_num, dimension) - 1));
		}
	}
	return weight;
}

double exact_independent_normalization_constant(NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, int dimension, bool posterior_bool) {
	double val = 0;
	for (int combination_number = 0; combination_number < 2; combination_number++) {
		double un_normalized_wight = un_normalized_independent_crf_weight(dimension, combination_number, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, posterior_bool);
		val += exp(un_normalized_wight);
	}
	return log(val);
}

double exact_independent_probability(double normalization_constant, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, int combination_number, int dimension, bool posterior_bool) {
	double prob = exp(un_normalized_independent_crf_weight(dimension, combination_number, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, posterior_bool) - normalization_constant);
	return prob;
}


double exact_independent_marginal_probability(double normalization_constant, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, int dimension, bool posterior_bool) {
	double marginal_prob = exact_independent_probability(normalization_constant, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, 1, dimension, posterior_bool);
	return marginal_prob;
}



// [[Rcpp::export]]
List update_independent_marginal_probabilities_exact_inference_cpp(NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int number_of_pairs, bool posterior_bool) {
	// Rcpp::Rcout << discrete_outliers(0,1) << std::endl;  
	// bool temp = discrete_outliers(0,1) == discrete_outliers(0,1);
	// Rcpp::Rcout << temp << std::endl;  

	// Initialize output matrices
	NumericMatrix probabilities(feat.nrow(), number_of_dimensions);
	NumericMatrix probabilities_pairwise(feat.nrow(), number_of_pairs);

	// Loop through samples
	for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
		// Loop through dimensions
		for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
			double normalization_constant = exact_independent_normalization_constant(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, dimension, posterior_bool);
			probabilities(sample_num, dimension) = exact_independent_marginal_probability(normalization_constant, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, dimension, posterior_bool);
		}
		int dimension_counter = 0;
		for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
			for (int dimension2=dimension; dimension2 < number_of_dimensions; dimension2++) {
				if (dimension != dimension2) {
					probabilities_pairwise(sample_num, dimension_counter) = probabilities(sample_num, dimension)*probabilities(sample_num, dimension2);
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
double compute_independent_crf_likelihood_exact_inference_cpp(NumericMatrix posterior, NumericMatrix posterior_pairwise, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, double lambda, double lambda_pair, double lambda_singleton) {
	// Rcpp::Rcout << feat(1,1) << std::endl;  
	// Initialize output likelihood
	double log_likelihood = 0;
	// Loop through samples
	for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
		// Loop through dimensions
		int dimension = 0;
		int dimension_counter = 0;
		for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
			double normalization_constant = exact_independent_normalization_constant(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, dimension, false);
			log_likelihood = log_likelihood - normalization_constant;
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

double logistic_regression_normalization_constant(NumericMatrix feat, double intercept, NumericVector theta, int sample_num) {
	double val = 0;
	for (int combination_number = 0; combination_number < 2; combination_number++) {
		double un_normalized_weight = 0;
		for (int d = 0; d < feat.ncol(); d++) {
			un_normalized_weight += feat(sample_num, d)*theta(d)*combination_number;
		}
		un_normalized_weight += intercept*combination_number;
		val += exp(un_normalized_weight);
	}
	return log(val);
}


// [[Rcpp::export]]
double compute_logistic_regression_likelihood_exact_inference_cpp(NumericMatrix posterior, NumericMatrix feat, double intercept, NumericVector theta, double lambda) {
	// Rcpp::Rcout << feat(1,1) << std::endl;  
	// Initialize output likelihood
	double log_likelihood = 0;
	// Loop through samples
	for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
		double normalization_constant = logistic_regression_normalization_constant(feat, intercept, theta, sample_num);
		log_likelihood = log_likelihood - normalization_constant;
		for (int d = 0; d < feat.ncol(); d++) {
			log_likelihood += theta(d)*feat(sample_num, d)*posterior(sample_num, 0);
		}
		log_likelihood += intercept*posterior(sample_num, 0);
	}

	log_likelihood = log_likelihood/feat.nrow();

	// Add L2 penalties
	for (int d = 0; d < feat.ncol(); d++) {
		log_likelihood = log_likelihood - .5*lambda*(theta(d)*theta(d));
	}
	return log_likelihood;
}

// [[Rcpp::export]]
NumericMatrix logistic_regression_predictions(NumericMatrix feat, double intercept, NumericVector theta) {
	NumericMatrix probabilities(feat.nrow(), 1);

	// Loop through samples
	double un_normalized_weight = 0;
	for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
		un_normalized_weight = 0;
		double normalization_constant = logistic_regression_normalization_constant(feat, intercept, theta, sample_num);
		for (int d = 0; d < feat.ncol(); d++) {
			un_normalized_weight += theta(d)*feat(sample_num, d);
		}
		un_normalized_weight += intercept;
		probabilities(sample_num,0) = exp(un_normalized_weight - normalization_constant);
	}
	return probabilities;
}

