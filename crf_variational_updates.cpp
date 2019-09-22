#include <vector>
#include <math.h>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <time.h>
#include <Rcpp.h>
using namespace Rcpp;





double variational_update(int sample_num, int dimension, NumericMatrix feat, NumericMatrix discrete_outliers, NumericMatrix probabilities, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, bool posterior_bool) {
	double term_b = 0.0;
	// Add singleton (intercept) term
	double term_a = theta_singleton(dimension);
	// Add pairwise term
	int dimension_counter = 0;
	for (int dimension1 = 0; dimension1 < number_of_dimensions; dimension1++) {
		for (int dimension2=dimension1; dimension2 < number_of_dimensions; dimension2++) {
			if (dimension1 != dimension2) {
				if (dimension1 == dimension) {
					term_a += theta_pair(0, dimension_counter)*probabilities(sample_num, dimension2);
				} else if (dimension2 == dimension) {
					term_a += theta_pair(0, dimension_counter)*probabilities(sample_num, dimension1);

				}
				dimension_counter += 1;
			}
		}
	}
	// Add feature based term
	for (int d = 0; d < feat.ncol(); d++) {
		term_a += feat(sample_num, d)*theta(d, dimension);
	}
	// Check to see if we are supposed to incorperate expression data && whether the expression data is observed
	if (posterior_bool == true && discrete_outliers(sample_num, dimension) == discrete_outliers(sample_num, dimension)) {
		term_a += log(phi_outlier(dimension, discrete_outliers(sample_num, dimension) - 1));
		term_b += log(phi_inlier(dimension, discrete_outliers(sample_num, dimension) - 1));
	}
	double variational_prob = exp(term_a)/(exp(term_a) + exp(term_b));
	return variational_prob;
}

double compute_elbo(NumericMatrix mu, int sample_num, int number_of_dimensions, NumericMatrix feat, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta) {
	// Loop through dimensions
	int dimension = 0;
	int dimension_counter = 0;
	double elbo = 0;
	double constant = 0;
	for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
		elbo += theta_singleton(dimension)*mu(sample_num, dimension);
		if (mu(sample_num, dimension) != 0.0 && mu(sample_num,dimension) != 1.0) {
			elbo -= (mu(sample_num, dimension))*log(constant + mu(sample_num, dimension));
			elbo -= (1.0 - mu(sample_num, dimension))*log(constant + 1.0 - mu(sample_num, dimension));
		}
		for (int d = 0; d < feat.ncol(); d++) {
			elbo += theta(d, dimension)*feat(sample_num, d)*mu(sample_num, dimension);
		}
		for (int dimension2=dimension; dimension2 < number_of_dimensions; dimension2++) {
			if (dimension != dimension2) {
				for (int theta_pair_dimension=0; theta_pair_dimension < theta_pair.nrow(); theta_pair_dimension++) {
					if (theta_pair_dimension == 0) {
						elbo += theta_pair(theta_pair_dimension, dimension_counter)*mu(sample_num, dimension)*mu(sample_num,dimension2);
					} else {
						elbo += feat(sample_num, (theta_pair_dimension-1))*theta_pair(theta_pair_dimension, dimension_counter)*mu(sample_num, dimension)*mu(sample_num,dimension2);
					}
				}
				dimension_counter += 1;
			}
		}
	}
	return elbo;

}

NumericMatrix variational_optimization(NumericMatrix probabilities, int sample_num, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, bool posterior_bool, double step_size, double convergence_thresh, std::mt19937 g) {
	// std::random_device rd;
	// std::mt19937 g(rd());
	double diff_prob = 0;
	int iteration_counter = 0;
	int convergence = 0;
	NumericMatrix prev_prob(1, number_of_dimensions);
	srand(0);
	while (convergence == 0) {
		diff_prob = 0;
		std::vector<int> v(number_of_dimensions);
		for (int dimension=0; dimension < number_of_dimensions; dimension++) {
			v[dimension] = dimension;
		}
 		std::shuffle(v.begin(), v.end(), g);

		for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
			// Rcpp:Rcout << v[dimension] << std::endl;
			// Save previous probability in order to check for convergence
			prev_prob(0, v[dimension]) = probabilities(sample_num, v[dimension]);
			// Update probability according to coordinate descent updates
			probabilities(sample_num, v[dimension]) = (1.0-step_size)*probabilities(sample_num, v[dimension]) + step_size*variational_update(sample_num, v[dimension], feat, discrete_outliers, probabilities, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, posterior_bool);
			// Compute differences
			diff_prob += std::abs(prev_prob(0, v[dimension]) - probabilities(sample_num, v[dimension]));
		}

		if (diff_prob/((double)number_of_dimensions) < convergence_thresh) {
			convergence = 1;
		}
		if (iteration_counter > 15000) {
			Rcpp::Rcout << "SKIPPED" << std::endl;  
			convergence = 1;
		}
		iteration_counter += 1;
	}
	// Rcpp::Rcout << iteration_counter << std::endl;  
	return probabilities;

}

double logistic_regression_initialization(int sample_num, int dimension, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, bool posterior_bool) {
	double term_b = 0.0;
	// Add singleton (intercept) term
	double term_a = theta_singleton(dimension);
	// Add pairwise term
	int dimension_counter = 0;

	// Add feature based term
	for (int d = 0; d < feat.ncol(); d++) {
		term_a += feat(sample_num, d)*theta(d, dimension);
	}
	// Check to see if we are supposed to incorperate expression data && whether the expression data is observed
	if (posterior_bool == true && discrete_outliers(sample_num, dimension) == discrete_outliers(sample_num, dimension)) {
		term_a += log(phi_outlier(dimension, discrete_outliers(sample_num, dimension) - 1));
		term_b += log(phi_inlier(dimension, discrete_outliers(sample_num, dimension) - 1));
	}
	double variational_prob = exp(term_a)/(exp(term_a) + exp(term_b));
	return variational_prob;



}

// [[Rcpp::export]]
List update_marginal_probabilities_vi_cpp(NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int number_of_pairs, double step_size, double convergence_thresh, NumericMatrix probability_init, bool posterior_bool) {
	std::random_device rd;
    std::mt19937 g(rd());
	//Rcpp::Rcout << discrete_outliers(0,1) << std::endl;  
	// bool temp = discrete_outliers(0,1) == discrete_outliers(0,1);
	// Rcpp::Rcout << temp << std::endl;  

	// Initialize output matrices
	NumericMatrix probabilities(feat.nrow(), number_of_dimensions);
	//NumericMatrix probabilities_rand(feat.nrow(), number_of_dimensions);
	//NumericMatrix probabilities_prev(feat.nrow(), number_of_dimensions);
	// NumericMatrix probabilities_log_reg_init(feat.nrow(), number_of_dimensions);
	NumericMatrix probabilities_pairwise(feat.nrow(), number_of_pairs);

	// Initialize temp matrix
	NumericMatrix prev_prob(1, number_of_dimensions);
	double diff_prob = 0;
	int convergence = 0;
	int iteration_counter = 0;
	//double elbo_start = 0;
	//double elbo = 0;
	//double elbo_rand = 0;
	//double elbo_prev = 0;
	//double elbo_total = 0;
	// double elbo_log_reg_init = 0;
	//double hit_counter = 0;
	// Loop through samples
	for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
		//Rcpp::Rcout << sample_num << std::endl;  

		// Initialize posterior prob
		for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
			//probabilities_prev(sample_num, dimension) = probability_init(sample_num, dimension);
			probabilities(sample_num, dimension) = ((double) rand() / (RAND_MAX));
			// probabilities_log_reg_init(sample_num, dimension) = logistic_regression_initialization(sample_num, dimension, feat, discrete_outliers, theta_singleton, theta, phi_inlier, phi_outlier, posterior_bool);
		}



		convergence = 0;
		iteration_counter = 0;


		probabilities = variational_optimization(probabilities, sample_num, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, posterior_bool, step_size, convergence_thresh, g);
		// probabilities = compute_elbo(probabilities, sample_num, number_of_dimensions, feat, theta_singleton, theta_pair, theta);

		//probabilities_prev = variational_optimization(probabilities_prev, sample_num, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, posterior_bool, step_size, convergence_thresh, g);
		//elbo_prev = compute_elbo(probabilities_prev, sample_num, number_of_dimensions, feat, theta_singleton, theta_pair, theta);


		// Compute pairwise probabilities (just product of marginals here)
		int dimension_counter = 0;
		for (int dimension1 = 0; dimension1 < number_of_dimensions; dimension1++) {
			for (int dimension2=dimension1; dimension2 < number_of_dimensions; dimension2++) {
				if (dimension1 != dimension2) {
					probabilities_pairwise(sample_num, dimension_counter) = probabilities(sample_num, dimension1)*probabilities(sample_num, dimension2);
					dimension_counter += 1;
				}
			}
		}
		// elbo = compute_elbo(probabilities, sample_num, number_of_dimensions, feat, theta_singleton, theta_pair, theta);
		// elbo_total = elbo_total + elbo;

	}
	//Rcpp::Rcout << elbo_total << std::endl;  
	List ret;
	ret["probability"] = probabilities;
	ret["probability_pairwise"] = probabilities_pairwise;
	return ret;
}


// [[Rcpp::export]]
double compute_crf_likelihood_vi_cpp(NumericMatrix posterior, NumericMatrix posterior_pairwise, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, double lambda, double lambda_pair, double lambda_singleton, NumericMatrix mu, NumericMatrix mu_pairwise) {
	double constant = 0;
	// Rcpp::Rcout << feat(1,1) << std::endl;  
	// Initialize output likelihood
	double log_likelihood = 0;
	// Loop through samples
	for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
		// Loop through dimensions
		int dimension = 0;
		int dimension_counter = 0;
		for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
			log_likelihood += theta_singleton(dimension)*posterior(sample_num, dimension) - theta_singleton(dimension)*mu(sample_num, dimension);
			if (mu(sample_num, dimension) != 0.0 && mu(sample_num,dimension) != 1.0) {
				log_likelihood += (mu(sample_num, dimension))*log(constant + mu(sample_num, dimension));
				log_likelihood += (1.0 - mu(sample_num, dimension))*log(constant + 1.0 - mu(sample_num, dimension));
			}
			for (int d = 0; d < feat.ncol(); d++) {
				log_likelihood += theta(d, dimension)*feat(sample_num, d)*posterior(sample_num, dimension) - theta(d, dimension)*feat(sample_num, d)*mu(sample_num, dimension);
			}
			for (int dimension2=dimension; dimension2 < number_of_dimensions; dimension2++) {
				if (dimension != dimension2) {
					for (int theta_pair_dimension=0; theta_pair_dimension < theta_pair.nrow(); theta_pair_dimension++) {
						if (theta_pair_dimension == 0) {
							log_likelihood += theta_pair(theta_pair_dimension, dimension_counter)*posterior_pairwise(sample_num, dimension_counter) - theta_pair(theta_pair_dimension, dimension_counter)*mu_pairwise(sample_num, dimension_counter);
						} else {
							log_likelihood += feat(sample_num, (theta_pair_dimension-1))*theta_pair(theta_pair_dimension, dimension_counter)*posterior_pairwise(sample_num, dimension_counter) - feat(sample_num, (theta_pair_dimension-1))*theta_pair(theta_pair_dimension, dimension_counter)*mu_pairwise(sample_num, dimension_counter);
						}
					}
					dimension_counter += 1;
				}
			}
		}
	}

	log_likelihood = (log_likelihood/feat.nrow());

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


