#include <vector>
#include <math.h>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <time.h>
#include <Rcpp.h>
#include <random>
using namespace Rcpp;

// Update variational probability according to coordinate descent updates for one (sample, dimension) pair
double variational_update(int sample_num, int dimension, NumericMatrix feat, NumericMatrix discrete_outliers, NumericMatrix probabilities, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, bool posterior_bool) {
  // Iniitialize term_a and term_b
  // term_b corresponds to the un-normalized weight of the mean field approximation to the current dimension assuming it takes on a value of 0
  // term_a corresponds to the un-normalized weight of the mean field approximation to the current dimension assuming it takes on a value of 1
  double term_b = 0.0;
  // Add singleton (intercept) term
  double term_a = theta_singleton(dimension);
  // Add pairwise term for each dimension that shares an edge with the current dimension
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
  // Normalize the terms to get a probability
  double variational_prob = exp(term_a)/(exp(term_a) + exp(term_b));
  return variational_prob;
}

// Use Mean Field Variational inference to infer posterior probabilities according to Conditional Random Field for this sample
NumericMatrix variational_optimization(NumericMatrix probabilities, int sample_num, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, bool posterior_bool, double step_size, double convergence_thresh, std::mt19937 g) {
  // Initialize temporary variables
  double diff_prob = 0;  // Keeps track of average absolute difference in mean field estimates between this iteration and the previous iteration
  int iteration_counter = 0; // Keeps track of the number of iterations
  int convergence = 0;  // 0 if not converged, 1 if converged
  NumericMatrix prev_prob(1, number_of_dimensions);  //Vector keep track of last iterations mean field estimates
  
  // Iterate until converged
  while (convergence == 0) {
    // Initialize variable the keeps track of average difference in mean field estimates between this iteration and the previous iteration to zero
    diff_prob = 0;
    
    // Create vector v that shuffles the dimensions so we can loop through dimensions in random ordering
    std::vector<int> v(number_of_dimensions);
    for (int dimension=0; dimension < number_of_dimensions; dimension++) {
      v[dimension] = dimension;
    }
    std::shuffle(v.begin(), v.end(), g);
    
    // Iterate through dimensions
    for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
      // Save previous probability in order to check for convergence
      prev_prob(0, v[dimension]) = probabilities(sample_num, v[dimension]);
      // Update variational probability according to coordinate descent updates for one (sample, dimension) pair
      probabilities(sample_num, v[dimension]) = (1.0-step_size)*probabilities(sample_num, v[dimension]) + step_size*variational_update(sample_num, v[dimension], feat, discrete_outliers, probabilities, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, posterior_bool);
      // Compute absolute difference between this iteration and the previous iteration
      diff_prob += std::abs(prev_prob(0, v[dimension]) - probabilities(sample_num, v[dimension]));
    }
    
    // Check for convergence via convergence threshold
    if (diff_prob/((double)number_of_dimensions) < convergence_thresh) {
      convergence = 1;
    }
    // Check for convergence via max iteration threshold
    if (iteration_counter > 15000) {
      Rcpp::Rcout << "SKIPPED" << std::endl;  
      convergence = 1;
    }
    
    iteration_counter += 1;
  }
  return probabilities;
  
}

// Compute both marginal Posterior probability as well as marginal pairwise probabilities for Conditional random field dimensions using Mean Field Variational Inference
// If posterior_bool==true, compute P(Z|E,G)
// If posterior_bool==false, compute P(Z|G)
// [[Rcpp::export]]
List update_marginal_probabilities_vi_cpp(NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int number_of_pairs, double step_size, double convergence_thresh, NumericMatrix probability_init, bool posterior_bool) {
  // Initialize random number generator
  std::random_device rd;
  std::mt19937 g(rd());
  
  // Initialize output matrices to keep tract of variational posteriors and variational pairwise posteriors
  NumericMatrix probabilities(feat.nrow(), number_of_dimensions);
  NumericMatrix probabilities_pairwise(feat.nrow(), number_of_pairs);
  
  // Iterate through samples
  for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
    
    // Initialize variational probabilities in each dimension to a random number between 0 and 1
    for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
      probabilities(sample_num, dimension) = ((double) rand() / (RAND_MAX));
    }
    
    // Use Mean Field Variational inference to infer posterior probabilities according to Conditional Random Field for this sample
    probabilities = variational_optimization(probabilities, sample_num, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, posterior_bool, step_size, convergence_thresh, g);
    
    // Compute pairwise probabilities (just product of marginals here as we are using mean field and all marginals are independent)
    int dimension_counter = 0;
    for (int dimension1 = 0; dimension1 < number_of_dimensions; dimension1++) {
      for (int dimension2=dimension1; dimension2 < number_of_dimensions; dimension2++) {
        if (dimension1 != dimension2) {
          probabilities_pairwise(sample_num, dimension_counter) = probabilities(sample_num, dimension1)*probabilities(sample_num, dimension2);
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
