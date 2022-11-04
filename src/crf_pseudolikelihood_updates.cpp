#include <vector>
#include <math.h>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <time.h>
#include <Rcpp.h>
using namespace Rcpp;

// Calculate the un-normalized pseudolikelihood weight for this value
double un_normalized_pseudolikelihood_crf_weight(int dimension, int combination_number, NumericMatrix feat, NumericMatrix posterior, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, bool posterior_bool) {
  // Initialize weight
  double weight = 0;
  // Add contribution from intercept
  weight += combination_number*theta_singleton(dimension);
  // Add contribution from feature weights
  for (int d = 0; d < feat.ncol(); d++) {
    weight += combination_number*feat(sample_num,d)*theta(d,dimension);
  }
  // Loop through pairs of dimensions
  int dimension_counter = 0;
  for (int dimension1 = 0; dimension1 < number_of_dimensions; dimension1++) {
    for (int dimension2=dimension1; dimension2 < number_of_dimensions; dimension2++) {
      if (dimension1 != dimension2) {
        // Add contribution from edge weights
        if (dimension1 == dimension) {
          weight += theta_pair(0, dimension_counter)*posterior(sample_num, dimension2)*combination_number;
        } else if (dimension2 == dimension) {
          weight += theta_pair(0, dimension_counter)*posterior(sample_num, dimension1)*combination_number;
        }
        dimension_counter += 1;
      }
    }
  }
  // Check to see if we are supposed to incorperate expression data && whether the expression data is observed
  if (posterior_bool == true && discrete_outliers(sample_num, dimension) == discrete_outliers(sample_num, dimension)) {
    if (combination_number == 1) {
      weight += log(phi_outlier(dimension, discrete_outliers(sample_num, dimension) - 1));
    } else {
      weight += log(phi_inlier(dimension, discrete_outliers(sample_num, dimension) - 1));
    }
  }
  return weight;
}

// Extract normalization constant for this (sample, dimension) pair
double exact_pseudolikelihood_normalization_constant(NumericMatrix feat, NumericMatrix discrete_outliers, NumericMatrix posterior, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, int dimension, bool posterior_bool) {
  // Initialize normalization constant
  double val = 0;
  // Loop through all possible values this dimension can take on according to pseudolikelihood (ie. 2 values: 0 or 1)
  for (int combination_number = 0; combination_number < 2; combination_number++) {
    // Calculate the un-normalized pseudolikelihood weight for this value
    double un_normalized_wight = un_normalized_pseudolikelihood_crf_weight(dimension, combination_number, feat, posterior, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, posterior_bool);
    val += exp(un_normalized_wight);
  }
  return log(val);
}

// Compute probabilitiy P(Z_dimension=1) for the current sample according to the pseudolikelihood
// If posterior_bool == true, compute P(Z_dimension = 1 | G,E)
// If posterior_bool == false, compute P(Z_dimension = 1 | G)
double exact_pseudolikelihood_marginal_probability(double normalization_constant, NumericMatrix feat, NumericMatrix discrete_outliers, NumericMatrix posterior, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int sample_num, int dimension, bool posterior_bool) {	
  int combination_number = 1;
  double marginal_prob = exp(un_normalized_pseudolikelihood_crf_weight(dimension, combination_number, feat, posterior, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, posterior_bool) - normalization_constant);
  return marginal_prob;
}


// Compute both marginal Posterior probability as well as marginal pairwise probabilities for Pseudolikelihood-based Conditional random field using exact inference
// If posterior_bool==true, compute P(Z|E,G)
// If posterior_bool==false, compute P(Z|G)
// [[Rcpp::export]]
List update_pseudolikelihood_marginal_probabilities_exact_inference_cpp(NumericMatrix feat, NumericMatrix discrete_outliers, NumericMatrix posterior, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int number_of_pairs, bool posterior_bool) {
  // Initialize output matrices
  NumericMatrix probabilities(feat.nrow(), number_of_dimensions);
  // Note: There are two pairwise probabilities as we are using psuedolikelihood
  NumericMatrix probabilities_pairwise1(feat.nrow(), number_of_pairs);
  NumericMatrix probabilities_pairwise2(feat.nrow(), number_of_pairs);
  
  // Loop through samples
  for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
    // Loop through dimensions
    for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
      // Extract normalization constant for this (sample, dimension) pair
      double normalization_constant = exact_pseudolikelihood_normalization_constant(feat, discrete_outliers, posterior, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, dimension, posterior_bool);
      // Compute probabilitiy P(Z_dimension=1) for the current sample according to the pseudolikelihood
      probabilities(sample_num, dimension) = exact_pseudolikelihood_marginal_probability(normalization_constant, feat, discrete_outliers, posterior, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, dimension, posterior_bool);
    }
    // Iterate through pairs of dimensions
    int dimension_counter = 0;
    for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
      for (int dimension2=dimension; dimension2 < number_of_dimensions; dimension2++) {
        if (dimension != dimension2) {
          // For this dimension pair, compute pairwise probabilities according to pseudolikelihood
          probabilities_pairwise1(sample_num, dimension_counter) = probabilities(sample_num, dimension)*posterior(sample_num, dimension2);
          probabilities_pairwise2(sample_num, dimension_counter) = posterior(sample_num, dimension)*probabilities(sample_num, dimension2);
          dimension_counter += 1;
        }
      }
      
    }
  }
  
  List ret;
  ret["probability"] = probabilities;
  ret["probability_pairwise1"] = probabilities_pairwise1;
  ret["probability_pairwise2"] = probabilities_pairwise2;
  return ret;
}

// Compute exact pseudo-likelihood of K=number_of_dimensions dimensionsal Conditional Random Field (CRF)
// [[Rcpp::export]]
double compute_pseudolikelihood_crf_likelihood_exact_inference_cpp(NumericMatrix posterior, NumericMatrix posterior_pairwise, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, double lambda, double lambda_pair, double lambda_singleton) {
  // Initialize output likelihood
  double log_likelihood = 0;
  // Loop through samples
  for (int sample_num = 0; sample_num < feat.nrow(); sample_num++) {
    // Loop through dimensions
    //int dimension_counter = 0;
    for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
      // Compute normalization constant for a particular (sample, dimension) pair
      double normalization_constant = exact_pseudolikelihood_normalization_constant(feat, discrete_outliers, posterior, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, sample_num, dimension, false);
      // Subtract normalization constant from log likelihood
      log_likelihood = log_likelihood - normalization_constant;
      // Add intercept to log likelihood
      log_likelihood += theta_singleton(dimension)*posterior(sample_num, dimension);
      // Add feature terms to genomic annotations
      for (int d = 0; d < feat.ncol(); d++) {
        log_likelihood += theta(d, dimension)*feat(sample_num, d)*posterior(sample_num, dimension);
      }
      // Loop through pairs of dimensions
      int dimension_counter = 0;
      for (int dimension1 = 0; dimension1 < number_of_dimensions; dimension1++) {
        for (int dimension2=dimension1; dimension2 < number_of_dimensions; dimension2++) {
          // Add edge weights for a particular dimension pair (edge)
          if (dimension1 != dimension2) {
            if (dimension1 == dimension) {
              log_likelihood += theta_pair(0, dimension_counter)*posterior(sample_num, dimension2)*posterior(sample_num, dimension1);
            } else if (dimension2 == dimension) {
              log_likelihood += theta_pair(0, dimension_counter)*posterior(sample_num, dimension1)*posterior(sample_num, dimension2);
            }
            dimension_counter += 1;
          }
        }
      }
    }
  }
  
  // Normalize log likelihood by the number of sample
  log_likelihood = log_likelihood/feat.nrow();
  
  // Add L2 penalties
  int dimension_counter = 0;
  for (int dimension = 0; dimension < number_of_dimensions; dimension++) {
    // L2 prior on intercepts. NOTE: Generally lambda_singleton is set to zero
    log_likelihood = log_likelihood - .5*lambda_singleton*(theta_singleton(dimension)*theta_singleton(dimension));
    // L2 prior on edge weights
    for (int dimension2=dimension; dimension2 < number_of_dimensions; dimension2++) {
      if (dimension != dimension2) {
        log_likelihood = log_likelihood - lambda_pair*(theta_pair(0,dimension_counter)*theta_pair(0, dimension_counter));
        dimension_counter += 1;
      }
    }
    // L2 prior n feature weights
    for (int d = 0; d < feat.ncol(); d++) {
      log_likelihood = log_likelihood - .5*lambda*(theta(d,dimension)*theta(d,dimension));
    }
  }
  return log_likelihood;
}
