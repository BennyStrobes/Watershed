// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// extract_all_binary_combinations
NumericMatrix extract_all_binary_combinations(int n);
RcppExport SEXP _WatershedR_extract_all_binary_combinations(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(extract_all_binary_combinations(n));
    return rcpp_result_gen;
END_RCPP
}
// update_marginal_probabilities_exact_inference_cpp
List update_marginal_probabilities_exact_inference_cpp(NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int number_of_pairs, bool posterior_bool);
RcppExport SEXP _WatershedR_update_marginal_probabilities_exact_inference_cpp(SEXP featSEXP, SEXP discrete_outliersSEXP, SEXP theta_singletonSEXP, SEXP theta_pairSEXP, SEXP thetaSEXP, SEXP phi_inlierSEXP, SEXP phi_outlierSEXP, SEXP number_of_dimensionsSEXP, SEXP number_of_pairsSEXP, SEXP posterior_boolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type feat(featSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type discrete_outliers(discrete_outliersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta_singleton(theta_singletonSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta_pair(theta_pairSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi_inlier(phi_inlierSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi_outlier(phi_outlierSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_dimensions(number_of_dimensionsSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_pairs(number_of_pairsSEXP);
    Rcpp::traits::input_parameter< bool >::type posterior_bool(posterior_boolSEXP);
    rcpp_result_gen = Rcpp::wrap(update_marginal_probabilities_exact_inference_cpp(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, number_of_pairs, posterior_bool));
    return rcpp_result_gen;
END_RCPP
}
// compute_crf_likelihood_exact_inference_cpp
double compute_crf_likelihood_exact_inference_cpp(NumericMatrix posterior, NumericMatrix posterior_pairwise, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, double lambda, double lambda_pair, double lambda_singleton);
RcppExport SEXP _WatershedR_compute_crf_likelihood_exact_inference_cpp(SEXP posteriorSEXP, SEXP posterior_pairwiseSEXP, SEXP featSEXP, SEXP discrete_outliersSEXP, SEXP theta_singletonSEXP, SEXP theta_pairSEXP, SEXP thetaSEXP, SEXP phi_inlierSEXP, SEXP phi_outlierSEXP, SEXP number_of_dimensionsSEXP, SEXP lambdaSEXP, SEXP lambda_pairSEXP, SEXP lambda_singletonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type posterior(posteriorSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type posterior_pairwise(posterior_pairwiseSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type feat(featSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type discrete_outliers(discrete_outliersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta_singleton(theta_singletonSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta_pair(theta_pairSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi_inlier(phi_inlierSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi_outlier(phi_outlierSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_dimensions(number_of_dimensionsSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_pair(lambda_pairSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_singleton(lambda_singletonSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_crf_likelihood_exact_inference_cpp(posterior, posterior_pairwise, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, lambda, lambda_pair, lambda_singleton));
    return rcpp_result_gen;
END_RCPP
}
// update_pseudolikelihood_marginal_probabilities_exact_inference_cpp
List update_pseudolikelihood_marginal_probabilities_exact_inference_cpp(NumericMatrix feat, NumericMatrix discrete_outliers, NumericMatrix posterior, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int number_of_pairs, bool posterior_bool);
RcppExport SEXP _WatershedR_update_pseudolikelihood_marginal_probabilities_exact_inference_cpp(SEXP featSEXP, SEXP discrete_outliersSEXP, SEXP posteriorSEXP, SEXP theta_singletonSEXP, SEXP theta_pairSEXP, SEXP thetaSEXP, SEXP phi_inlierSEXP, SEXP phi_outlierSEXP, SEXP number_of_dimensionsSEXP, SEXP number_of_pairsSEXP, SEXP posterior_boolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type feat(featSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type discrete_outliers(discrete_outliersSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type posterior(posteriorSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta_singleton(theta_singletonSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta_pair(theta_pairSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi_inlier(phi_inlierSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi_outlier(phi_outlierSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_dimensions(number_of_dimensionsSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_pairs(number_of_pairsSEXP);
    Rcpp::traits::input_parameter< bool >::type posterior_bool(posterior_boolSEXP);
    rcpp_result_gen = Rcpp::wrap(update_pseudolikelihood_marginal_probabilities_exact_inference_cpp(feat, discrete_outliers, posterior, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, number_of_pairs, posterior_bool));
    return rcpp_result_gen;
END_RCPP
}
// compute_pseudolikelihood_crf_likelihood_exact_inference_cpp
double compute_pseudolikelihood_crf_likelihood_exact_inference_cpp(NumericMatrix posterior, NumericMatrix posterior_pairwise, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, double lambda, double lambda_pair, double lambda_singleton);
RcppExport SEXP _WatershedR_compute_pseudolikelihood_crf_likelihood_exact_inference_cpp(SEXP posteriorSEXP, SEXP posterior_pairwiseSEXP, SEXP featSEXP, SEXP discrete_outliersSEXP, SEXP theta_singletonSEXP, SEXP theta_pairSEXP, SEXP thetaSEXP, SEXP phi_inlierSEXP, SEXP phi_outlierSEXP, SEXP number_of_dimensionsSEXP, SEXP lambdaSEXP, SEXP lambda_pairSEXP, SEXP lambda_singletonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type posterior(posteriorSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type posterior_pairwise(posterior_pairwiseSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type feat(featSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type discrete_outliers(discrete_outliersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta_singleton(theta_singletonSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta_pair(theta_pairSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi_inlier(phi_inlierSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi_outlier(phi_outlierSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_dimensions(number_of_dimensionsSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_pair(lambda_pairSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_singleton(lambda_singletonSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_pseudolikelihood_crf_likelihood_exact_inference_cpp(posterior, posterior_pairwise, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, lambda, lambda_pair, lambda_singleton));
    return rcpp_result_gen;
END_RCPP
}
// update_marginal_probabilities_vi_cpp
List update_marginal_probabilities_vi_cpp(NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int number_of_pairs, double step_size, double convergence_thresh, NumericMatrix probability_init, bool posterior_bool);
RcppExport SEXP _WatershedR_update_marginal_probabilities_vi_cpp(SEXP featSEXP, SEXP discrete_outliersSEXP, SEXP theta_singletonSEXP, SEXP theta_pairSEXP, SEXP thetaSEXP, SEXP phi_inlierSEXP, SEXP phi_outlierSEXP, SEXP number_of_dimensionsSEXP, SEXP number_of_pairsSEXP, SEXP step_sizeSEXP, SEXP convergence_threshSEXP, SEXP probability_initSEXP, SEXP posterior_boolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type feat(featSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type discrete_outliers(discrete_outliersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta_singleton(theta_singletonSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta_pair(theta_pairSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi_inlier(phi_inlierSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi_outlier(phi_outlierSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_dimensions(number_of_dimensionsSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_pairs(number_of_pairsSEXP);
    Rcpp::traits::input_parameter< double >::type step_size(step_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type convergence_thresh(convergence_threshSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type probability_init(probability_initSEXP);
    Rcpp::traits::input_parameter< bool >::type posterior_bool(posterior_boolSEXP);
    rcpp_result_gen = Rcpp::wrap(update_marginal_probabilities_vi_cpp(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, number_of_pairs, step_size, convergence_thresh, probability_init, posterior_bool));
    return rcpp_result_gen;
END_RCPP
}
// update_independent_marginal_probabilities_exact_inference_cpp
List update_independent_marginal_probabilities_exact_inference_cpp(NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, int number_of_pairs, bool posterior_bool);
RcppExport SEXP _WatershedR_update_independent_marginal_probabilities_exact_inference_cpp(SEXP featSEXP, SEXP discrete_outliersSEXP, SEXP theta_singletonSEXP, SEXP theta_pairSEXP, SEXP thetaSEXP, SEXP phi_inlierSEXP, SEXP phi_outlierSEXP, SEXP number_of_dimensionsSEXP, SEXP number_of_pairsSEXP, SEXP posterior_boolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type feat(featSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type discrete_outliers(discrete_outliersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta_singleton(theta_singletonSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta_pair(theta_pairSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi_inlier(phi_inlierSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi_outlier(phi_outlierSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_dimensions(number_of_dimensionsSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_pairs(number_of_pairsSEXP);
    Rcpp::traits::input_parameter< bool >::type posterior_bool(posterior_boolSEXP);
    rcpp_result_gen = Rcpp::wrap(update_independent_marginal_probabilities_exact_inference_cpp(feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, number_of_pairs, posterior_bool));
    return rcpp_result_gen;
END_RCPP
}
// compute_independent_crf_likelihood_exact_inference_cpp
double compute_independent_crf_likelihood_exact_inference_cpp(NumericMatrix posterior, NumericMatrix posterior_pairwise, NumericMatrix feat, NumericMatrix discrete_outliers, NumericVector theta_singleton, NumericMatrix theta_pair, NumericMatrix theta, NumericMatrix phi_inlier, NumericMatrix phi_outlier, int number_of_dimensions, double lambda, double lambda_pair, double lambda_singleton);
RcppExport SEXP _WatershedR_compute_independent_crf_likelihood_exact_inference_cpp(SEXP posteriorSEXP, SEXP posterior_pairwiseSEXP, SEXP featSEXP, SEXP discrete_outliersSEXP, SEXP theta_singletonSEXP, SEXP theta_pairSEXP, SEXP thetaSEXP, SEXP phi_inlierSEXP, SEXP phi_outlierSEXP, SEXP number_of_dimensionsSEXP, SEXP lambdaSEXP, SEXP lambda_pairSEXP, SEXP lambda_singletonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type posterior(posteriorSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type posterior_pairwise(posterior_pairwiseSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type feat(featSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type discrete_outliers(discrete_outliersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta_singleton(theta_singletonSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta_pair(theta_pairSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi_inlier(phi_inlierSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi_outlier(phi_outlierSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_dimensions(number_of_dimensionsSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_pair(lambda_pairSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_singleton(lambda_singletonSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_independent_crf_likelihood_exact_inference_cpp(posterior, posterior_pairwise, feat, discrete_outliers, theta_singleton, theta_pair, theta, phi_inlier, phi_outlier, number_of_dimensions, lambda, lambda_pair, lambda_singleton));
    return rcpp_result_gen;
END_RCPP
}
// compute_logistic_regression_likelihood_exact_inference_cpp
double compute_logistic_regression_likelihood_exact_inference_cpp(NumericMatrix posterior, NumericMatrix feat, double intercept, NumericVector theta, double lambda);
RcppExport SEXP _WatershedR_compute_logistic_regression_likelihood_exact_inference_cpp(SEXP posteriorSEXP, SEXP featSEXP, SEXP interceptSEXP, SEXP thetaSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type posterior(posteriorSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type feat(featSEXP);
    Rcpp::traits::input_parameter< double >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_logistic_regression_likelihood_exact_inference_cpp(posterior, feat, intercept, theta, lambda));
    return rcpp_result_gen;
END_RCPP
}
// logistic_regression_predictions
NumericMatrix logistic_regression_predictions(NumericMatrix feat, double intercept, NumericVector theta);
RcppExport SEXP _WatershedR_logistic_regression_predictions(SEXP featSEXP, SEXP interceptSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type feat(featSEXP);
    Rcpp::traits::input_parameter< double >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(logistic_regression_predictions(feat, intercept, theta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_WatershedR_extract_all_binary_combinations", (DL_FUNC) &_WatershedR_extract_all_binary_combinations, 1},
    {"_WatershedR_update_marginal_probabilities_exact_inference_cpp", (DL_FUNC) &_WatershedR_update_marginal_probabilities_exact_inference_cpp, 10},
    {"_WatershedR_compute_crf_likelihood_exact_inference_cpp", (DL_FUNC) &_WatershedR_compute_crf_likelihood_exact_inference_cpp, 13},
    {"_WatershedR_update_pseudolikelihood_marginal_probabilities_exact_inference_cpp", (DL_FUNC) &_WatershedR_update_pseudolikelihood_marginal_probabilities_exact_inference_cpp, 11},
    {"_WatershedR_compute_pseudolikelihood_crf_likelihood_exact_inference_cpp", (DL_FUNC) &_WatershedR_compute_pseudolikelihood_crf_likelihood_exact_inference_cpp, 13},
    {"_WatershedR_update_marginal_probabilities_vi_cpp", (DL_FUNC) &_WatershedR_update_marginal_probabilities_vi_cpp, 13},
    {"_WatershedR_update_independent_marginal_probabilities_exact_inference_cpp", (DL_FUNC) &_WatershedR_update_independent_marginal_probabilities_exact_inference_cpp, 10},
    {"_WatershedR_compute_independent_crf_likelihood_exact_inference_cpp", (DL_FUNC) &_WatershedR_compute_independent_crf_likelihood_exact_inference_cpp, 13},
    {"_WatershedR_compute_logistic_regression_likelihood_exact_inference_cpp", (DL_FUNC) &_WatershedR_compute_logistic_regression_likelihood_exact_inference_cpp, 5},
    {"_WatershedR_logistic_regression_predictions", (DL_FUNC) &_WatershedR_logistic_regression_predictions, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_WatershedR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
