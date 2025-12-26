#' @title Lasso-Based Prognostic Index Calculation
#' @name Lasso_PI_scores_f
#'
#' @details The function first checks the validity of input data. It renames
#' the survival time and event columns in the training and test datasets,
#' then fits a Cox proportional hazards model using Lasso regularization
#' via cross-validation. Based on the Lasso model, the function calculates
#' the prognostic index (PI) for both the training and test datasets,
#' appending the PI values to the original data frames.
#' 
#' @param train_data A data frame containing the training data, including
#' survival time, event status, and features.
#' @param test_data A data frame containing the test data, including survival
#' time, event status, and features.
#' @param nfolds An integer specifying the number of cross-validation folds
#' for LASSO.
#' @param col_num An integer indicating the starting column index for the
#' feature variables.
#' @param surv_time A character string specifying the column name for
#' survival time in the data.
#' @param surv_event A character string specifying the column name for
#' survival event status in the data.
#' @param n_repeats An integer specifying how many times repeated CV should be run to assess feature selection stability.
#' @param alpha A numeric value (0–1) specifying the mixing parameter for Elastic Net (1 = LASSO, 0 < alpha < 1 = Elastic Net).
#' @param freq_threshold Numeric value between 0 and 1. 
#'   Minimum selection frequency required to retain a feature across
#'   repeated LASSO runs (e.g., 0.6 keeps features selected in ≥60% of repeats).
#' @return A list with the following components:
#' \itemize{
#'   \item \code{Train_Lasso_key_variables}: A data frame of selected key
#'   variables and their coefficients.
#'   \item \code{Train_PI_data}: A data frame with training data, including
#'   survival time, event status, and PI.
#'   \item \code{Test_PI_data}: A data frame with test data, including
#'   survival time, event status, and PI.
#'   \item \code{cvfit}: A `cv.glmnet` object representing the fitted LASSO
#'   Cox model.
#'   \item \code{Stability_table}: A data frame showing each feature and the
#'   proportion of repeated CV runs in which it was selected (selection frequency).
#' }
#'
#' @import survival
#' @import survminer
#' @import ggplot2
#' @import glmnet
#' @importFrom grDevices gray
#' @importFrom stats as.formula coef na.omit predict
#'
#' @examples
#' data(Train_Norm_data, package = "CPSM")
#' data(Test_Norm_data, package = "CPSM")
#' Result_PI <- Lasso_PI_scores_f(
#' train_data = Train_Norm_data,
#' test_data = Test_Norm_data,
#' nfolds = 5,
#' col_num = 21,
#' n_repeats = 10,       
#' alpha = 1, 
#' freq_threshold = 0.6,      
#' surv_time = "OS_month",
#' surv_event = "OS")
#'
#' @export





Lasso_PI_scores_f <- function(train_data, test_data, nfolds, col_num,
                              surv_time, surv_event, n_repeats,  alpha, freq_threshold = 0.6) {
  # Check if any input variable is empty
  if (length(train_data) == 0 || length(test_data) == 0 ||
    length(nfolds) == 0 || length(col_num) == 0 ||
    length(surv_time) == 0 || length(surv_event) == 0) {
    message("Error: Empty input variable detected.")
  }

  # Check if nfolds is valid
  if (is.na(nfolds) || nfolds <= 0) {
    message("Error: nfolds is missing or invalid.")
  }
  # Check if col_num is valid
  if (is.na(col_num) || col_num <= 0) {
    message("Error: col_num is missing or invalid.")
  }
  if (!surv_time %in% colnames(train_data)) {
    message("Error: surv_time column not found in the data.")
  }
  if (!surv_event %in% colnames(train_data)) {
    message("Error: surv_event column not found in the data.")
  }


  # load data
  tr_data1 <- train_data
  te_data1 <- test_data

  # rename survival time and event column name
  colnames(tr_data1)[colnames(tr_data1) == surv_time] <- "OS_month"
  colnames(tr_data1)[colnames(tr_data1) == surv_event] <- "OS"

  # rename survival time and event column name
  colnames(te_data1)[colnames(te_data1) == surv_time] <- "OS_month"
  colnames(te_data1)[colnames(te_data1) == surv_event] <- "OS"

  # LASSO COX
  # create survival object
  surv_object <- survival::Surv(time = tr_data1$OS_month, event = tr_data1$OS)


# Repeated CV for feature stability
# -------------------------
	feature_names <- colnames(tr_data1)[col_num:ncol(tr_data1)]
	selection_matrix <- matrix(
  	0,
  	nrow = length(feature_names),
  	ncol = n_repeats,
  	dimnames = list(feature_names, paste0("Rep_", 1:n_repeats))
	)

	for (r in 1:n_repeats) {

  	cvfit_tmp <- glmnet::cv.glmnet(
    	as.matrix(tr_data1[col_num:ncol(tr_data1)]),
   	 surv_object,
    	family = "cox",
    	type.measure = "C",
    	nfolds = min(nfolds, floor(nrow(tr_data1)/2)), # safer for small n
    	alpha = alpha,
   	 maxit = 1000
  	)

  	coef_tmp <- coef(cvfit_tmp, s = cvfit_tmp$lambda.min)
  	coef_vec <- as.numeric(coef_tmp)

  	selection_matrix[, r] <- as.integer(coef_vec != 0)
	}

	# Compute selection frequency
	selection_frequency <- rowMeans(selection_matrix)
	Stability_table <- data.frame(
  	Feature = names(selection_frequency),
  	Selection_Frequency = round(selection_frequency, 3)
	)

	# Keep only features selected in >=60% of repeats
	stable_features <- Stability_table$Feature[
  	Stability_table$Selection_Frequency >= freq_threshold
	]

	if (length(stable_features) == 0) {
    	warning("No stable features selected; using all features instead.")
    	stable_features <- colnames(tr_data1)[col_num:ncol(tr_data1)]
	}
#-----------------------------------------------#

  # develop model to select key features using LASSO method and user defined
  # cross-folds
  cvfit1 <- glmnet::cv.glmnet(
    #as.matrix(tr_data1[col_num:ncol(tr_data1)]),
    as.matrix(tr_data1[, stable_features, drop = FALSE]),
    surv_object, # create survival object from the data
    family = "cox", # specify Cox PH model
    type.measure = "C",
    nfolds = nfolds,
    alpha = alpha,
    maxit = 1000
  )

  # calculate lambda minimum
  est.coef <- coef(cvfit1, s = cvfit1$lambda.min)
  # returns the p length coefficient vector
  # of the solution corresponding to lambda
  est.coef1 <- as.numeric(est.coef)
  active.k <- which(est.coef1 != 0)

  # Extract coefficinet values
  active.k.vals <- est.coef[active.k]

  key_variables <- as.data.frame(est.coef[est.coef[, 1] != 0, ])
  colnames(key_variables) <- c("coeff")
  key_variables <- round(key_variables, 3)

  # Create PI Index for training data
  sel_features2 <- key_variables
  sel_train2 <- as.data.frame(tr_data1[, colnames(tr_data1) %in%
    c(row.names(sel_features2)), ])

  # Make final files with selected features  and survival information
  train_feature_mat2 <- cbind(tr_data1["OS_month"], tr_data1["OS"], sel_train2)

  # Create prognostic index (PI)
  tr_PI <- train_feature_mat2[3:ncol(train_feature_mat2)]
  E <- length(tr_PI)

  PI_tr <- 0
  for (i in seq(from = 1, to = E, by = 1)) {
    PI_tr <- PI_tr + ((tr_PI[, i]) * (sel_features2[i, 1]))
  }

  # add PI as new column to the data
  tr_PI$PI <- PI_tr
  # combines PI information with survival info
  train_PI <- cbind(
    train_feature_mat2["OS"],
    train_feature_mat2["OS_month"], tr_PI
  )
  # PI for Test data
  sel_test <- as.data.frame(te_data1[, colnames(te_data1) %in%
    c(row.names(sel_features2)), ])

  # Make final files with selected features & survival info
  test_feature_mat2 <- cbind(te_data1["OS_month"], te_data1["OS"], sel_test)

  # Create prognostic index
  te_PI <- test_feature_mat2[3:ncol(test_feature_mat2)]

  E <- length(te_PI)
  sel_features2[2, 1]

  PI_te <- 0
  for (i in seq(from = 1, to = E, by = 1)) {
    PI_te <- PI_te + ((te_PI[, i]) * (sel_features2[i, 1]))
  }

  # add PI as new column to the data
  te_PI$PI <- PI_te
  test_PI <- as.data.frame(te_PI$PI)
  rownames(test_PI) <- rownames(te_PI)
  colnames(test_PI) <- c("PI")

  test_PI <- cbind(
    test_feature_mat2["OS"], test_feature_mat2["OS_month"],
    te_PI
  )
  # save selected data with PI value
  Train_Lasso_key_variables <- key_variables
  Train_PI_data <- train_PI
  Test_PI_data <- test_PI

  # Return a list containing data.
  return(list(
    Train_Lasso_key_variables = Train_Lasso_key_variables,
    Train_PI_data = Train_PI_data,
    Test_PI_data = Test_PI_data,
    Stability_table = Stability_table, 
    cvfit = cvfit1 # Return the cvfit object to plot externally
  ))
}
