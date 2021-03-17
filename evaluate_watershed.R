library(optparse)
library(PRROC)
source("watershed.R")

# Helper function to remove NAs from vector
remove_na <- function(x) {
	return(x[!is.na(x)])
}

#######################################
# Extract precision recall curves for test set (in each dimension seperately) using:
#### 1. Watershed predictions
#### 2. GAM predictions
#######################################
compute_auc_across_dimensions <- function(number_of_dimensions, watershed_posteriors, gam_posteriors, real_valued_outliers_test1, binary_outliers_test2) {
	# Initialize objects to store auc information
	auc_object_across_dimensions <- list()
	pos_list <- c()
	neg_list <- c()
  	# Loop through dimensions
  	for (dimension in 1:number_of_dimensions) {
  		# Pseudo gold standard
  		test_outlier_status <- binary_outliers_test2[,dimension]

		# Watershed evaluation curves
  		watershed_roc_obj <- roc.curve(scores.class0 = remove_na(watershed_posteriors[,dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]), scores.class1 = remove_na(watershed_posteriors[,dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]), curve = T)
  		watershed_pr_obj <- pr.curve(scores.class0 = remove_na(watershed_posteriors[,dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]), scores.class1 = remove_na(watershed_posteriors[,dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]), curve = T)
  	
  		# GAM evaluation curves
  		gam_roc_obj <- roc.curve(scores.class0 = remove_na(gam_posteriors[,dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]), scores.class1 = remove_na(gam_posteriors[,dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]), curve = T)
   		gam_pr_obj <- pr.curve(scores.class0 = remove_na(gam_posteriors[,dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[,dimension])]), scores.class1 = remove_na(gam_posteriors[,dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[,dimension])]), curve = T)

		evaROC <-	
		 list(watershed_sens=watershed_roc_obj$curve[,2],
              watershed_spec=1-watershed_roc_obj$curve[,1],
         	  watershed_auc=watershed_roc_obj$auc,
         	  watershed_pr_auc=watershed_pr_obj$auc.integral,
         	  watershed_recall=watershed_pr_obj$curve[,1],
         	  watershed_precision=watershed_pr_obj$curve[,2],
         	  GAM_sens=gam_roc_obj$curve[,2],
              GAM_spec=1-gam_roc_obj$curve[,1],
              GAM_auc=gam_roc_obj$auc,
         	  GAM_pr_auc=gam_pr_obj$auc.integral,
         	  GAM_recall=gam_pr_obj$curve[,1],
         	  GAM_precision=gam_pr_obj$curve[,2])

		auc_object_across_dimensions[[dimension]] <- list(evaROC=evaROC)
	}

	return(auc_object_across_dimensions)
}



#######################################
## Train model on non-N2 pairs and evaluate model on N2-pairs
#######################################
evaluate_watershed_shell <- function(input_file, number_of_dimensions, model_name, pseudoc, lambda_init, output_stem, n2_pair_pvalue_fraction, binary_pvalue_threshold, lambda_costs, nfolds, vi_step_size, vi_threshold) {
	########################
	## Load in data
	########################
	data_input <- load_watershed_data(input_file, number_of_dimensions, n2_pair_pvalue_fraction, binary_pvalue_threshold)
	# Parse data_input for evaluation-related relevent fields
	feat_all <- data_input$feat
	discrete_outliers_all <- data_input$outliers_discrete
	binary_outliers_all <- data_input$outliers_binary
	fraction_binary_outliers_all <- data_input$fraction_outliers_binary
	N2_pairs <- data_input$N2_pairs

	## Extract training data (ie. non-N2-pairs)
	feat_train <- feat_all[is.na(N2_pairs),]
	discrete_outliers_train <- as.matrix(discrete_outliers_all[is.na(N2_pairs),])
	binary_outliers_train <-  as.matrix(binary_outliers_all[is.na(N2_pairs),])

	## Extract test data (ie N2 pairs)
	#(has to be done seperately for RIVER vs Watershed)
	if (number_of_dimensions == 1) {
		# Extraction of test data for RIVER
		feat_test <- rbind(feat_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], feat_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
  		discrete_outliers_test1 <- as.matrix(c(discrete_outliers_all[!is.na(N2_pairs)][seq(from=1,to=sum(!is.na(N2_pairs)),by=2)], discrete_outliers_all[!is.na(N2_pairs)][seq(from=2,to=sum(!is.na(N2_pairs)),by=2)]))
  		discrete_outliers_test2 <- as.matrix(c(discrete_outliers_all[!is.na(N2_pairs)][seq(from=2,to=sum(!is.na(N2_pairs)),by=2)], discrete_outliers_all[!is.na(N2_pairs)][seq(from=1,to=sum(!is.na(N2_pairs)),by=2)]))
		binary_outliers_test1 <- as.matrix(c(binary_outliers_all[!is.na(N2_pairs)][seq(from=1,to=sum(!is.na(N2_pairs)),by=2)], binary_outliers_all[!is.na(N2_pairs)][seq(from=2,to=sum(!is.na(N2_pairs)),by=2)]))
  		binary_outliers_test2 <- as.matrix(c(fraction_binary_outliers_all[!is.na(N2_pairs)][seq(from=2,to=sum(!is.na(N2_pairs)),by=2)], fraction_binary_outliers_all[!is.na(N2_pairs)][seq(from=1,to=sum(!is.na(N2_pairs)),by=2)]))
		# Absolute pvalues from test prediction data set (to be used for RNA-only analysis)
  		real_valued_outliers_test1 <- -log10(abs(as.matrix(c(data_input$outlier_pvalues[!is.na(N2_pairs)][seq(from=1,to=sum(!is.na(N2_pairs)),by=2)], data_input$outlier_pvalues[!is.na(N2_pairs)][seq(from=2,to=sum(!is.na(N2_pairs)),by=2)]))) + 1e-7)

	} else {
		# Extraction of test data for Watershed
  		feat_test <- rbind(feat_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], feat_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
  		discrete_outliers_test1 <- rbind(discrete_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], discrete_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
  		discrete_outliers_test2 <- rbind(discrete_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),], discrete_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),])
  		binary_outliers_test1 <- rbind(binary_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], binary_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
  		binary_outliers_test2 <- rbind(fraction_binary_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),], fraction_binary_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),])
  		# Absolute pvalues from test prediction data set (to be used for RNA-only analysis)
  		real_valued_outliers_test1 <- -log10(abs(rbind(data_input$outlier_pvalues[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], data_input$outlier_pvalues[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])) + 1e-7)
  	}

	#######################################
	## Standardize Genomic Annotations (features)
	#######################################
	mean_feat <- apply(feat_all, 2, mean)
	sd_feat <- apply(feat_all, 2, sd)
 	feat_all <- scale(feat_all, center=mean_feat, scale=sd_feat)
 	feat_train <- scale(feat_train, center=mean_feat, scale=sd_feat)
 	feat_test <- scale(feat_test, center=mean_feat, scale=sd_feat)

 	#######################################
	## Fit Genomic Annotation Model (GAM)
	#######################################
	gam_data <- logistic_regression_genomic_annotation_model_cv(feat_train, binary_outliers_train, nfolds, lambda_costs, lambda_init)
	# Report optimal lambda learned from cross-validation data (if applicable)
	if (is.na(lambda_init)) {
		cat(paste0(nfolds,"-fold cross validation on GAM yielded optimal lambda of ", gam_data$lambda, "\n"))
	}
	# Compute GAM Predictions on test data in CPP file ("independent_crf_exact_updates.cpp")
	gam_posterior_test_obj <- update_independent_marginal_probabilities_exact_inference_cpp(feat_test, binary_outliers_test1, gam_data$gam_parameters$theta_singleton, gam_data$gam_parameters$theta_pair, gam_data$gam_parameters$theta, matrix(0,2,2), matrix(0,2,2), number_of_dimensions, choose(number_of_dimensions, 2), FALSE)
	gam_test_posteriors <- gam_posterior_test_obj$probability

	#######################################
	### Initialize phi using GAM
	#######################################
	# Compute GAM Predictions on training data in CPP file ("independent_crf_exact_updates.cpp")
	gam_posterior_train_obj <- update_independent_marginal_probabilities_exact_inference_cpp(feat_train, binary_outliers_test1, gam_data$gam_parameters$theta_singleton, gam_data$gam_parameters$theta_pair, gam_data$gam_parameters$theta, matrix(0,2,2), matrix(0,2,2), number_of_dimensions, choose(number_of_dimensions, 2), FALSE)
	gam_train_posteriors <- gam_posterior_train_obj$probability
	# Initialize Phi using GAM posteriors
	# ie. Compute MAP estimates of the coefficients defined by P(outlier_status| FR)
	phi_init <- map_phi_initialization(discrete_outliers_train, gam_train_posteriors, number_of_dimensions, pseudoc)

	#######################################
	### Fit Watershed Model
	#######################################
	watershed_model <- train_watershed_model(feat_train, discrete_outliers_train, phi_init, gam_data$gam_parameters$theta_pair, gam_data$gam_parameters$theta_singleton, gam_data$gam_parameters$theta, pseudoc, gam_data$lambda, number_of_dimensions, model_name, vi_step_size, vi_threshold)

	#######################################
	## Compute Watershed Posterior probabilities for held-out test data (ie the N2 pairs)
	#######################################
	posterior_info_test <- update_marginal_posterior_probabilities(feat_test, discrete_outliers_test1, watershed_model)
	posterior_prob_test <- posterior_info_test$probability  # Marginal posteriors
	posterior_pairwise_prob_test <- posterior_info_test$probability_pairwise  # Pairwise posteriors

	#######################################
	# Extract precision recall curves for test set (in each dimension seperately) using:
	#### 1. Watershed predictions
	#### 2. GAM predictions
	#######################################
	auc_object_across_dimensions <- compute_auc_across_dimensions(number_of_dimensions, posterior_prob_test, gam_test_posteriors, real_valued_outliers_test1, binary_outliers_test2)

	return(list(auc=auc_object_across_dimensions, model_params=watershed_model, gam_model_params=gam_data))
}







#########################################
# Command line arguments
#########################################
arguments <- parse_args(OptionParser(usage = "%prog [options]", description="Watershed command line args",
	option_list=list(
		make_option(c("-i","--input"), default = NULL, help="The Watershed input file  [default %default]"),
		make_option(c("-d","--number_dimensions"), default=1, help="The number of outlier variables."),
		make_option(c("-m","--model_name"), default="Watershed_exact", help="Name of model. Options are Watershed_exact, Watershed_approximate, and RIVER"),
		make_option(c("-p","--dirichlet_prior_parameter"), default=10, help="Parameter defining Dirichlet distribution that acts as a prior a Phi (the model parameters defining E|Z"),
		make_option(c("-l","--l2_prior_parameter"), default=.01, help="Parameter defining L2 (gaussian) distribution that acts as a prior on the parameters of the conditional random Field (the model parameters defining Z|G"),
		make_option(c("-o","--output_prefix"), default="watershed", help="Prefix of file name to save results to"),
		make_option(c("-n","--n2_pair_pvalue_fraction"), default=.1, help="Fraction of outliers (based on rank) that are considered an outlier for N2 pair analysis (this is done so each outlier type/signal has approximately the same distribution of positive outlier examples)"),
		make_option(c("-b","--binary_pvalue_threshold"), default=.1, help="Absolute p-value threshold used to create binary outliers used for Genomic Annotation Model"))
))
# process command line args
input_file <- arguments$input
number_of_dimensions <- arguments$number_dimensions
model_name <- arguments$model_name
pseudoc <- arguments$dirichlet_prior_parameter
lambda_init <- arguments$l2_prior_parameter
if (lambda_init == "NA") {
	lambda_init <- NA
}
output_stem <- arguments$output_prefix
n2_pair_pvalue_fraction <- arguments$n2_pair_pvalue_fraction
binary_pvalue_threshold <- arguments$binary_pvalue_threshold

# Change model to RIVER if there is only 1 dimension.
if (number_of_dimensions==1 & model_name != "RIVER") {
	cat("Only RIVER can be run on data with 1 dimension.\n Changing model to RIVER.\n")
	model_name <- "RIVER"
}

#########################################
# Fixed Parameters
#########################################
# If arguments$l2_prior_parameter == NA, perform grid search over the following values of lambda to determine optimal lambda
lambda_costs <- c(.1, .01, 1e-3)
# If arguments$l2_prior_parameter == NA, Number of folds to be used in K-fold cross validation for Genomic annotation model ()
nfolds <- 5
# Parameters used for Variational Optimization (only applies if arguments$model_name=="Watershed_approximate")
vi_step_size <- .8
vi_threshold <- 1e-8

# Set seed for reproducability (feel free to change!)
set.seed(1)


#######################################
## Train model on non-N2 pairs and evaluate model on N2-pairs
#######################################
evaluation_object <- evaluate_watershed_shell(input_file, number_of_dimensions, model_name, pseudoc, lambda_init, output_stem, n2_pair_pvalue_fraction, binary_pvalue_threshold, lambda_costs, nfolds, vi_step_size, vi_threshold)


#######################################
## Save evaluation object as .rds file
#######################################
saveRDS(evaluation_object, paste0(output_stem, "_evaluation_object.rds"))


#######################################
## Print area under precision-recall curves in each dimension
#######################################
for (dimension in 1:number_of_dimensions) {
	model_pr_auc <- evaluation_object$auc[[dimension]]$evaROC$watershed_pr_auc
	gam_pr_auc <- evaluation_object$auc[[dimension]]$evaROC$GAM_pr_auc
	cat(paste0(model_name, " area under precision-recall curve in dimension ", dimension,": ", model_pr_auc, "\n"))
	cat(paste0("GAM area under precision-recall curve in dimension ", dimension,": ", gam_pr_auc, "\n"))
}