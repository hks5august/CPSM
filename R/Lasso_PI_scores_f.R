#' This function will create PI (Prognostic Index) score based on selected
#' LASSO genes for training and test data.  Here, firstly, LASSO will select
#' genes with beta coefficient values based on COX lasso regression using
#' 5-fold cross validation. Subsequently, PI score will be calculated by
#' multiplying expression of genes with their beta coeff values.
#' @param train_data:args1 - training data (Patients data with clinical and
#' gene expression, where samples are in rows and features/genes are in
#' columns)
#' @param test_data:args2 - training data (Patients data with clinical and
#' gene expression, where samples are in rows and features/genes are in
#' columns)
#' @param nfolds:args3 - Number of folds for cross in cvglmnet model to select
#' top features
#' @param col_num:args4 - column number in data at where clinical info ends
#' @param surv_time:arg5 - name of column which contain survival time (in days)
#' information
#' @param surv_event:arg6 - name of column which contain survival
#' eventinformation
#' @import dplyr
#' @import survival
#' @import survminer
#' @import ggplot2
#' @import glmnet
#' @import svglite
#' @examples
#' data(Train_Norm_data, package = "CPSM")
#' data(Test_Norm_data, package = "CPSM")
#' Lasso_PI_scores_f(
#'   train_data = Train_Norm_data, test_data = Test_Norm_data,
#'   nfolds = 5, col_num = 21, surv_time = "OS_month", surv_event = "OS"
#' )
#' Usage:Lasso_PI_scores_f(
#'   train_data, test_data, nfolds, col_num, surv_time,
#'   surv_event
#' )
#' @export

Lasso_PI_scores_f <- function(train_data, test_data, nfolds, col_num,
                              surv_time, surv_event) {
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

  # develop model to select key features using LASSO method and user defined
  # cross-folds
  cvfit1 <- glmnet::cv.glmnet(
    as.matrix(tr_data1[col_num:ncol(tr_data1)]),
    surv_object, # create survival object from the data
    family = "cox", # specify Cox PH model
    type.measure = "C",
    nfolds = nfolds,
    alpha = 1,
    maxit = 1000
  )

  # calculate lambda minimum
  lambda_min <- cvfit1$lambda.min
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
    cvfit = cvfit1 # Return the cvfit object to plot externally
  ))
}
