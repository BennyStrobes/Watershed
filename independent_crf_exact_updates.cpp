#include <vector>
#include <math.h>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <time.h>
#include <Rcpp.h>
using namespace Rcpp;

// For a logistic regression outcome (combination number), compute the relative logistic regressin weight
double un_normalized_independent_crf_weight(int dimension, int combination_number, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, bool posterior_bool) {
	// Initialize weight
	double weight = 0;
	int dimension_counter = 0;
	// Contribution from intercept
	weight += combination_number*theta_singleton(dimension);
	// Contribution from features
	for (int d = 0; d < feat.ncol(); d++) {
		weight += combination_number*feat(sample_num,d)*theta(d,dimension);
	}
	// If posterior_bool is true and outlier measurement is observed, add contribution of E
	if (posterior_bool == true && discrete_outliers(sample_num, dimension) == discrete_outliers(sample_num, dimension)) {
		// Z == 1
		if (combination_number == 1) {
			weight += log(phi_outlier(dimension, discrete_outliers(sample_num, dimension) - 1));
		// Z == 0
		} else {
			weight += log(phi_inlier(dimension, discrete_outliers(sample_num, dimension) - 1));
		}
	}
	return weight;
}

// Compute logistic normalization constant
double exact_independent_normalization_constant(NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, int dimension, bool posterior_bool) {
	// Initialize variable to keep track of the normalization constant
	double val = 0;
	// Loop through all possible outcomes of logistic regression (only two outcomes)
	for (int combination_number = 0; combination_number < 2; combination_number++) {
		// For each outcome, calculate the relative weight
		double un_normalized_wight = un_normalized_independent_crf_weight(dimension, combination_number, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, posterior_bool);
		val += exp(un_normalized_wight);
	}
	return log(val);
}


// Compute logistic regression probability for a specific (sample, dimension) pair
double exact_independent_marginal_probability(double normalization_constant, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, int dimension, bool posterior_bool) {
	// Predict P(Z=1)
	int combination_number = 1;
	double marginal_prob = exp(un_normalized_independent_crf_weight(dimension, combination_number, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, posterior_bool) - normalization_constant);

	return marginal_prob;
}


// Compute both marginal Posterior probability latent variable as well as marginal pairwise probabilities
// If posterior_bool==true, compute P(Z|E,G)
// If posterior_bool==false, compute P(Z|G)
// [[Rcpp::export]]
List update_independent_marginal_probabilities_exact_inference_cpp(NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int number_of_pairs, bool posterior_bool) {
	// Initialize output matrices
	NumericMatrix probabilities(feat.nrow(), number_of_dimensions);
	NumericMatrix probabilities_pairwise(feat.nrow(), number_of_pairs);

	// Loop through samples
	for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
		// Loop through dimensions
		for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
			// Compute logistic normalization constant
			double normalization_constant = exact_independent_normalization_constant(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, dimension, posterior_bool);
			// Use normalization constant to compute marginal posterior probability 
			probabilities(sample_num, dimension) = exact_independent_marginal_probability(normalization_constant, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, dimension, posterior_bool);
		}
		// Loop through pairs of dimensions to compute pairwise probabilities
		// Because this is for k=number_of_dimensions independent logistic regression models, instead of an actual CRF, each probability is independent
		// And pairwise probability is just the product of the marginals
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

// Compute likelihood for K=number_of_dimensions independent logistic regression models
// [[Rcpp::export]]
double compute_independent_crf_likelihood_exact_inference_cpp(NumericMatrix posterior, NumericMatrix posterior_pairwise, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, double lambda, double lambda_pair, double lambda_singleton) {
	// Initialize output likelihood
	double log_likelihood = 0;
	// Loop through samples
	for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
		// Loop through dimensions
		int dimension = 0;
		int dimension_counter = 0;
		for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
			// Compute normalization constant for this (sample, dimension pair) and substract from log_likelihood
			double normalization_constant = exact_independent_normalization_constant(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, dimension, false);
			log_likelihood = log_likelihood - normalization_constant;
			// Add interecept (theta_singleton) to likelihood
			log_likelihood += theta_singleton(dimension)*posterior(sample_num, dimension);
			// Add contribution of features
			for (int d = 0; d < feat.ncol(); d++) {
				log_likelihood += theta(d, dimension)*feat(sample_num, d)*posterior(sample_num, dimension);
			}
		}
	}
	// Normalize likelihood by number of samples
	log_likelihood = log_likelihood/feat.nrow();

	// Add L2 penalties to likelihood
	int dimension_counter = 0;
	for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
		// Generally lambda_singleton==0, so this term does not contribute. But theoretically could
		log_likelihood = log_likelihood - .5*lambda_singleton*(theta_singleton(dimension)*theta_singleton(dimension));
		// Regularize feature weights
		for (int d = 0; d < feat.ncol(); d++) {
			log_likelihood = log_likelihood - .5*lambda*(theta(d,dimension)*theta(d,dimension));
		}
	}
	return log_likelihood;
}

// Calculate logistic regression normalization constant
double logistic_regression_normalization_constant(NumericMatrix feat, double intercept, NumericVector theta, int sample_num) {
	// Initialize variable to keep track of normalization constant
	double val = 0;
	// Loop through all possible outcomes of logistic regression (only two outcomes)
	for (int combination_number = 0; combination_number < 2; combination_number++) {
		// For each outcome, calculate relative contribution
		double un_normalized_weight = 0;
		// Contribution from features
		for (int d = 0; d < feat.ncol(); d++) {
			un_normalized_weight += feat(sample_num, d)*theta(d)*combination_number;
		}
		// Contribution from intercept
		un_normalized_weight += intercept*combination_number;
		// Add combination_number contribution to val
		val += exp(un_normalized_weight);
	}
	return log(val);
}


// Smaller version of "compute_independent_crf_likelihood_exact_inference_cpp" to just compute logistic regression likelihood for one dimensions
// [[Rcpp::export]]
double compute_logistic_regression_likelihood_exact_inference_cpp(NumericMatrix posterior, NumericMatrix feat, double intercept, NumericVector theta, double lambda) {
	// Rcpp::Rcout << feat(1,1) << std::endl;  
	// Initialize output likelihood
	double log_likelihood = 0;
	// Loop through samples
	for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
		// Calculate logistic regression normalization constant and substract from log likelihood
		double normalization_constant = logistic_regression_normalization_constant(feat, intercept, theta, sample_num);
		log_likelihood = log_likelihood - normalization_constant;
		// Add contribution of features to likelihood
		for (int d = 0; d < feat.ncol(); d++) {
			log_likelihood += theta(d)*feat(sample_num, d)*posterior(sample_num, 0);
		}
		// Add contribution of intercepts to likelihood
		log_likelihood += intercept*posterior(sample_num, 0);
	}

	// Normalize likelihood by the number of samples
	log_likelihood = log_likelihood/feat.nrow();

	// Add L2 penalties
	for (int d = 0; d < feat.ncol(); d++) {
		log_likelihood = log_likelihood - .5*lambda*(theta(d)*theta(d));
	}
	return log_likelihood;
}

// Make predictions P(Z=1) for logistic regression model
// [[Rcpp::export]]
NumericMatrix logistic_regression_predictions(NumericMatrix feat, double intercept, NumericVector theta) {
	// Initialize matrix keeping tracking of predictions across all samples
	NumericMatrix probabilities(feat.nrow(), 1);

	// Loop through samples
	double un_normalized_weight = 0;
	for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
		// Get normalization constant for this sample
		double normalization_constant = logistic_regression_normalization_constant(feat, intercept, theta, sample_num);
		// Get un-normalized weight for Z=1
		un_normalized_weight = 0;
		for (int d = 0; d < feat.ncol(); d++) {
			un_normalized_weight += theta(d)*feat(sample_num, d);
		}
		un_normalized_weight += intercept;
		// Normalize probability
		probabilities(sample_num,0) = exp(un_normalized_weight - normalization_constant);
	}
	return probabilities;
}

