library(lbfgs)
library(Rcpp)
sourceCpp("independent_crf_exact_updates.cpp")


# Convert outlier status into discretized random variables
get_discretized_outliers <- function(outlier_pvalues) {
	# initialize output
	outliers_discretized <- matrix(0,dim(outlier_pvalues)[1], dim(outlier_pvalues)[2])
	for (dimension in 1:ncol(outlier_pvalues)) {
		# Check if it is total expression
		if (min(outlier_pvalues[,dimension], na.rm=TRUE) < 0) {
			under_expression = outlier_pvalues[,dimension] < 0 & !is.na(outlier_pvalues[,dimension])
			log_pvalues = -log10(abs(outlier_pvalues[,dimension]) + 1e-6)
			log_pvalues[under_expression] = log_pvalues[under_expression]*-1
			#discretized <- cut(log_pvalues,breaks=c(-6.01,-4,-2,-1,1,2,4,6.01))
			discretized <- cut(log_pvalues, breaks=c(-6.01,-1,1,6.01))
		} else {
			log_pvalues = -log10(abs(outlier_pvalues[,dimension]) + 1e-6)
			# discretized <- cut(log_pvalues, 7)
			discretized <- cut(log_pvalues, breaks=c(-.01,1,4,6))
		}
		outliers_discretized[,dimension] = as.numeric(discretized)
	}
	colnames(outliers_discretized) = colnames(outlier_pvalues)
	return(outliers_discretized)
}


# Load in and parse Watershed input file
load_watershed_data <- function(input_file, number_of_dimensions, pvalue_fraction, pvalue_threshold) {
	raw_data <- read.table(input_file, header=TRUE)
	# Get genomic features (first 2 columns are line identifiers and last (number_of_dimensions+1) columns are outlier status' and N2 pair
	feat <- raw_data[,3:(ncol(raw_data)-number_of_dimensions-1)]
	# sample name as SubjectID:GeneName
	rownames(feat) <- paste(raw_data[,"SubjectID"], ":", raw_data[,"GeneName"],sep="")
	# Pvalues of outlier status of a particular sample (rows) for a particular outlier type (columns)
	outlier_pvalues <- as.matrix(raw_data[,(ncol(raw_data)-number_of_dimensions):(ncol(raw_data)-1)])
	# sample name as SubjectID:GeneName
	rownames(outlier_pvalues) <- paste(raw_data[,"SubjectID"], ":", raw_data[,"GeneName"],sep="")
	# Convert outlier status into binary random variables
	fraction_outliers_binary <- ifelse(abs(outlier_pvalues)<=.1,1,0) # Strictly for initialization of binary output matrix
	for (dimension_num in 1:number_of_dimensions) {
		ordered <- sort(abs(outlier_pvalues[,dimension_num]))
		max_val <- ordered[floor(length(ordered)*pvalue_fraction)]
		fraction_outliers_binary[,dimension_num] <- ifelse(abs(outlier_pvalues[,dimension_num])<=max_val,1,0)
	}
  	outliers_binary <- ifelse(abs(outlier_pvalues)<=pvalue_threshold,1,0)
	# Convert outlier status into discretized random variables
	outliers_discrete <- get_discretized_outliers(outlier_pvalues)
	# Extract array of N2 pairs
	N2_pairs=factor(raw_data[,"N2pair"], levels=unique(raw_data[,"N2pair"]))
	# Put all data into compact data structure
	data_input <- list(feat=as.matrix(feat), outlier_pvalues=as.matrix(outlier_pvalues),outliers_binary=as.matrix(outliers_binary), fraction_outliers_binary=as.matrix(fraction_outliers_binary),outliers_discrete=outliers_discrete, N2_pairs=N2_pairs)
	return(data_input)
}

# Compute log-likelihood of L2-regularized logistic regression model
compute_logistic_regression_likelihood <- function(x, y, feat, lambda) {
	intercept <- x[1]
	theta <- x[2:length(x)]
	# Compute likelihood in CPP file ("independent_crf_exact_updates.cpp")
	log_likelihood <- compute_logistic_regression_likelihood_exact_inference_cpp(y, feat, intercept, theta, lambda)
	return(-log_likelihood)
}

# Calculate gradient of L2-regularized logistic regression model
compute_logistic_regression_gradient <- function(x, y, feat, lambda) {
	# Extract intercept and coefficients
	intercept <- x[1]
	theta <- x[2:length(x)]

	# Make predictions in CPP file ("independent_crf_exact_updates.cpp")
	predictions <- logistic_regression_predictions(feat, intercept, theta)

	# Compute Gradient of singleton terms (intercepts)
	grad_singleton <- (colSums(y) - colSums(predictions))*(1/nrow(y))

	# Compute Gradient of theta terms (betas)
	grad_theta <- c()
	dimension <- 1
	temp_grad <- colSums(y[,dimension]*feat) - colSums(predictions[,dimension]*feat)
	grad_theta <- c(grad_theta, temp_grad)
	# add regularization
	grad_theta <- grad_theta*(1/nrow(y)) - lambda*theta

	grad <- c(grad_singleton, grad_theta)
	return(-grad)
}

# Compute number of edge pairs. In general is just N choose 2. But if N==1, we use the hack the number of edge_pairs is 1
get_number_of_edge_pairs <- function(number_of_dimensions) {
  val = 1
  if (number_of_dimensions > 1) {
    val = choose(number_of_dimensions,2)
  }
  return(val)
}


## Fit Genomic Annotation Model (GAM)
logistic_regression_genomic_annotation_model_cv <- function(feat_train, binary_outliers_train, nfolds, lambda_costs, lambda_init) {
	##################################
	# Some pre-processing
	##################################
	# Extract dimensionality of space from the data
	number_of_dimensions <- dim(binary_outliers_train)[2]
	number_of_features <- dim(feat_train)[2]

	# Initialize logistic regression model parameters to zeros
	gradient_variable_vec <- rep(0, number_of_features+1)
	
	#Randomly shuffle the data
	random_shuffling_indices <- sample(nrow(feat_train))
	feat_train_shuff <- as.matrix(feat_train[random_shuffling_indices,])
	binary_outliers_train_shuff <- as.matrix(binary_outliers_train[random_shuffling_indices,])

	##################################
	# Select value of lambda to use
	##################################
	# If lambda_init == NA, we will do K-fold cross validation to select the optimal lambda
	if (is.na(lambda_init)) {
		#Create nfolds equally size folds
		folds <- cut(seq(1,nrow(feat_train_shuff)),breaks=nfolds,labels=FALSE)
		# Initialize array to keep track of the average Area Under (Precision-Recall Curve) across different values of lambda (lambda_costs)
		avg_aucs <- c()

		# Iterate across lambdas in lambda_costs 
		for (cost_iter in 1:length(lambda_costs)) {
			lambda <- lambda_costs[cost_iter]
			# Initialize array to keep track of auc in each of the n-folds
			aucs <- c()
			# Loop through n-folds
			for(i in 1:nfolds){
				# Initialize array to keep track of logistic regression probabiliites in test samples (put in 'pos' array if we know test sample is positive via held out label. Other way around for 'neg' array)
				pos <- c()
				neg <- c()
    			#Segement your data into training and test for this fold
    			testIndexes <- which(folds==i,arr.ind=TRUE)
    			feat_test_fold <- feat_train_shuff[testIndexes,]
    			outliers_test_fold <- as.matrix(binary_outliers_train_shuff[testIndexes,])
    			feat_train_fold <- feat_train_shuff[-testIndexes,]
    			outliers_train_fold <- as.matrix(binary_outliers_train_shuff[-testIndexes,])

    			# Perform logistic regression in each dimension seperately
    			for (dimension in 1:number_of_dimensions) {
    				# Remove any samples with NA for this outlier dimension
    				observed_training_indices <- !is.na(outliers_train_fold[,dimension])
    				observed_training_outliers <- as.matrix(outliers_train_fold[observed_training_indices, dimension])
    				observed_training_feat <- feat_train_fold[observed_training_indices,]
    				observed_testing_indices <- !is.na(outliers_test_fold[,dimension])
    				observed_testing_outliers <- as.matrix(outliers_test_fold[observed_testing_indices, dimension])
    				observed_testing_feat <- feat_test_fold[observed_testing_indices,]
    			
    				# Optimize logistic gregression using LBFGS
    				lbfgs_output <- lbfgs(compute_logistic_regression_likelihood, compute_logistic_regression_gradient, gradient_variable_vec, y=observed_training_outliers, feat=observed_training_feat, lambda=lambda, invisible=1)
    			 
    				if (lbfgs_output$convergence != 0) {
    					print("ERRROR in logistic regression optimization!")
    					print(lbfgs_output$convergence)
    				}

    				# Make predictions on test data using learned Logistic regression model
    				predictions <- c(logistic_regression_predictions(observed_testing_feat, lbfgs_output$par[1], lbfgs_output$par[2:length(lbfgs_output$par)]))
    				# Add precictions to array
    				pos <- c(pos, predictions[observed_testing_outliers==1])
    				neg <- c(neg, predictions[observed_testing_outliers==0])
    			}

    			# Compute Precision recall curve for this fold
    			pr_obj <- pr.curve(scores.class0=pos, scores.class1=neg,curve=T)
    			# Get area under precision-recall curve
    			auc <- pr_obj$auc.integral
    			aucs <- c(aucs, auc)
			}
			# Compute the median across
			avg_aucs <- c(avg_aucs, median(aucs))
		}
		# Get best lambda (ie, the one with highest avg auc across folds)
		best_index <- which(avg_aucs==max(avg_aucs))[1]  # [1] for tie breakers
		best_lambda <- lambda_costs[best_index]
	# If lambda_init != NA, use the user-specified values
  	} else {
  		best_lambda = lambda_init
  	}

	##################################
	# Using optimal lambda, run GAM on full data to:
	## 1. Train GAM
	## 2. Use parameters in GAM to initialize watershed parameters
	##################################
	# Initialize parameter variables
	theta_pair = matrix(0 ,1, get_number_of_edge_pairs(number_of_dimensions))
	beta_init = matrix(0,number_of_features+1, number_of_dimensions)
	theta_singleton = beta_init[1,]
	theta = as.matrix(beta_init[2:(number_of_features + 1),])
	gam_parameters = list(theta_pair=theta_pair, theta_singleton=theta_singleton, theta=theta)

	# Run GAM in each dimension
  	for (dimension in 1:number_of_dimensions) {
  		# Remove any samples with NA for this outlier dimension
  		observed_training_indices <- !is.na(binary_outliers_train_shuff[,dimension])
  		observed_training_outliers <- as.matrix(binary_outliers_train_shuff[observed_training_indices, dimension])
  		observed_training_feat <- feat_train_shuff[observed_training_indices,]

  		# Run Logistic regression
  		lbfgs_output <- lbfgs(compute_logistic_regression_likelihood, compute_logistic_regression_gradient, gradient_variable_vec, y=observed_training_outliers, feat=observed_training_feat, lambda=best_lambda, invisible=1)

  		if (lbfgs_output$convergence != 0) {
    		print("ERRROR!")
    		print(lbfgs_output$convergence)
    	}
    	# Add optimal GAM parameters to data structure
    	gam_parameters$theta_singleton[dimension] <- lbfgs_output$par[1]
    	gam_parameters$theta[,dimension] <- lbfgs_output$par[2:length(lbfgs_output$par)]
  	}
  	return(list(lambda=best_lambda, gam_parameters=gam_parameters))
}


# Compute MAP estimates of the coefficients defined by P(outlier_status| FR)
map_phi_initialization <- function(discrete_outliers, posterior, number_of_dimensions, pseudoc) {
	num_bins = 3
	# Initialize output matrices
	phi_outlier <- matrix(1, number_of_dimensions, num_bins)  
	phi_inlier <- matrix(1, number_of_dimensions, num_bins)
	# Count number of times we fall into each bin
	for (bin_number in 1:num_bins) {
		phi_outlier[,bin_number] <- colSums(((discrete_outliers==bin_number)*posterior),na.rm=TRUE)
		phi_inlier[,bin_number] <- colSums(((discrete_outliers==bin_number)*(1-posterior)),na.rm=TRUE)
	}
	# Add prior
	for (dimension_number in 1:number_of_dimensions) {
		phi_outlier[dimension_number,] = phi_outlier[dimension_number,] + pseudoc
		phi_inlier[dimension_number,] = phi_inlier[dimension_number,] + pseudoc
	}
	# Normalize
	phi_outlier <- phi_outlier/rowSums(phi_outlier)
	phi_inlier <- phi_inlier/rowSums(phi_inlier)

	# Combine into compact object
	phi_init <- list(inlier_component = phi_inlier, outlier_component = phi_outlier)

	return(phi_init)
}


