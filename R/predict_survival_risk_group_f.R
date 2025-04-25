#' @title Survival Risk Group Prediction Function
#'
#' @description  This function builds a Random Forest (RF) model to classify samples into high-risk or low-risk survival groups
#' based on selected features.
#'
#' @name predict_survival_risk_group_f
#'
#' @param selected_train_data A data frame containing normalized expression data of selected features and survival information for
#' training set.
#' @param selected_test_data A data frame containing normalized expression data of selected features and survival information for
#' test set.

#' @param Feature_List A list of  features to be used for predicting the risk class in the
#' model.
#'
#' @return A list containing the following components:
#' \item{best_model}{The Best RF model.}
#' \item{Train_results}{Prediction resulst on Training data}
#' \item{Test_results}{Prediction resulst on Test data}
#' \item{misclassification_results}{Prediction results of best on model on training and Test data }
#' \item{All_tree_matrix}{Prediction results on training and test data based on all models using different number of trees }
#'
#' @importFrom caret confusionMatrix
#' @importFrom randomForestSRC rfsrc
#' @importFrom stats median complete.cases


#' @examples
#' # Example usage of the predict_survival_risk_group_f function
#' data(Train_PI_data, package = "CPSM")
#' data(Test_PI_data, package = "CPSM")
#' data(Key_PI_list , package = "CPSM")
#' Results_Risk_group_Prediction<-  predict_survival_risk_group_f(selected_train_data = Train_PI_data,
#' selected_test_data = Test_PI_data, Feature_List = Key_PI_list )
#'
#' @export



predict_survival_risk_group_f <-  function(selected_train_data, selected_test_data,
                                         Feature_List) {


  # Check if any input variable is empty
  if (length(selected_train_data) == 0 ||
      length( selected_test_data) == 0 ||
      length(Feature_List) == 0) {
    message("Error: Empty input variable detected.")
  }


   # Load input data
  Train <- selected_train_data
  Test <- selected_test_data


  # Check for required survival outcome columns
  required_columns <- c("OS", "OS_month")
  missing_train_cols <- setdiff(required_columns, colnames(Train))
  missing_test_cols <- setdiff(required_columns, colnames(Test))

  if (length(missing_train_cols) > 0 || length(missing_test_cols) > 0) {
    stop(paste("Missing required columns in training or test data:",
               paste(c(missing_train_cols, missing_test_cols), collapse = ", ")))
  }

  # Process selected features
  clinical_features <- Feature_List$ID  # Remove the first element
  clinical_features <- intersect(clinical_features, colnames(Train))
  clinical_features <- intersect(clinical_features, colnames(Test))

  #head(Train)
  #head(Test)

  # Convert character variables to factors
  for (feature in clinical_features) {
    if (is.character(Train[[feature]])) Train[[feature]] <- as.factor(Train[[feature]])
    if (is.character(Test[[feature]])) Test[[feature]] <- as.factor(Test[[feature]])
  }

  # Remove rows with missing values
  Train <- Train[complete.cases(Train[, clinical_features]), ]
  Test <- Test[complete.cases(Test[, clinical_features]), ]

  if (nrow(Train) == 0 | nrow(Test) == 0) {
    stop("Training or test data has no complete cases for selected features.")
  }

  # Define Actual Risk Groups based on median OS_month
  median_surv_train <- median(Train$OS_month, na.rm = TRUE)
  Train$Actual_Risk_Group <- factor(ifelse(Train$OS_month > median_surv_train, "Low_Risk", "High_Risk"))
  Test$Actual_Risk_Group <- factor(ifelse(Test$OS_month > median_surv_train, "Low_Risk", "High_Risk"))


  # ----- Classification Prediction model to Risk groups ------ #

  # Define ntree values
  ntree_values <- c(10, 20, 50, 100, 250, 500, 750, 1000)


  # Define a function to extract relevant metrics from confusion matrix
  extract_metrics <- function(cm) {
    if (!is.null(cm)) {
      TP <- cm$table[2, 2]
      TN <- cm$table[1, 1]
      FP <- cm$table[1, 2]
      FN <- cm$table[2, 1]

      accuracy <- round(cm$overall["Accuracy"] * 100, 2)
      sensitivity <- round(cm$byClass["Sensitivity"] * 100, 2)
      specificity <- round(cm$byClass["Specificity"] * 100, 2)
      precision <- round(TP / (TP + FP), 3)
      recall <- round(TP / (TP + FN), 3)
      F1 <- round(2 * (precision * recall) / (precision + recall), 3)
      FPR <- round(FP / (FP + TN), 3)
      FNR <- round(FN / (FN + TP), 3)

      return(c(accuracy, sensitivity, specificity, precision, recall, F1, FPR, FNR))
    } else {
      return(rep(NA, 8))
    }
  }

  # Initialize best model tracking
  best_ntree <- NULL
  best_test_accuracy <- 0
  best_model <- NULL

  # Initialize an empty matrix to store all metrics
  all_metrics_matrix <- data.frame()
  # Iterate over ntree values
  for (ntree in ntree_values) {
    cat("Training", "_with ntree =", ntree, "\n")

    # Define formula dynamically
    formula <- as.formula(paste("Actual_Risk_Group ~", paste(clinical_features, collapse = " + ")))

    # Train Random Forest model
    rf_model <- rfsrc(formula, data = Train, ntree = ntree, importance = TRUE)

    # Predictions
    train_pred <- predict(rf_model, newdata = Train)$class
    test_pred <- predict(rf_model, newdata = Test)$class

    # Compute confusion matrices
    #cm_train <- caret::compute_cm(train_pred, Train$Actual_Risk_Group)
    #cm_test <- caret::compute_cm(test_pred, Test$Actual_Risk_Group)


    # Compute confusion matrices
    cm_train <- confusionMatrix(as.factor(train_pred), as.factor(Train$Actual_Risk_Group))
    cm_test <- confusionMatrix(as.factor(test_pred), as.factor(Test$Actual_Risk_Group))


    # Extract metrics
    train_metrics <- extract_metrics(cm_train)
    test_metrics <- extract_metrics(cm_test)

    # Store all metrics in a single matrix
    metrics_matrix <- data.frame(
      ntree = ntree,
      Train_Accuracy = train_metrics[1],
      Train_Sensitivity = train_metrics[2],
      Train_Specificity = train_metrics[3],
      Train_Precision = train_metrics[4],
      Train_Recall = train_metrics[5],
      Train_F1 = train_metrics[6],
      Train_FPR = train_metrics[7],
      Train_FNR = train_metrics[8],
      Test_Accuracy = test_metrics[1],
      Test_Sensitivity = test_metrics[2],
      Test_Specificity = test_metrics[3],
      Test_Precision = test_metrics[4],
      Test_Recall = test_metrics[5],
      Test_F1 = test_metrics[6],
      Test_FPR = test_metrics[7],
      Test_FNR = test_metrics[8]
    )

    # Write to output file containing results for trees (ntree) for RF
    #write.table(metrics_matrix, file= output_file, sep = "\t", row.names = FALSE, col.names = !file.exists(output_file), append = TRUE)
    all_metrics_matrix <- rbind(all_metrics_matrix, metrics_matrix)


    # Save best model based on test accuracy
    if (!is.na(test_metrics[1]) && test_metrics[1] > best_test_accuracy) {
      best_ntree <- ntree
      best_test_accuracy <- test_metrics[1]
      best_model <- rf_model
    }
  }


  # Save predictions for the best model
  if (!is.null(best_model)) {

    # Create the file path
    #best_model_file <- file.path(Temp_path, paste0("Best_RF_Model_", feature_type, "_", best_ntree, ".rds"))

    # Save the best model
    #saveRDS(best_model, file = best_model_file)


    best_train_pred <- predict(best_model, newdata = Train)$class
    best_test_pred <- predict(best_model, newdata = Test)$class

    best_train_pred_prob <- predict(best_model, newdata = Train)$predicted
    best_test_pred_prob <- predict(best_model, newdata = Test)$predicted

    train_prob <- round(apply(best_train_pred_prob, 1, max), 3)
    test_prob <- round(apply(best_test_pred_prob, 1, max), 3)

    Train_results <- data.frame(Sample_ID = rownames(Train),  # Sample ID column
                                Actual = Train$Actual_Risk_Group,
                                Predicted_Risk_Group = best_train_pred,  # Add Predicted Risk Group
                                High_Risk_Prob = if ("High_Risk" %in% colnames(best_train_pred_prob)) best_train_pred_prob[, "High_Risk"] else rep(NA, nrow(best_train_pred_prob)),
                                Low_Risk_Prob = if ("Low_Risk" %in% colnames(best_train_pred_prob)) best_train_pred_prob[, "Low_Risk"] else rep(NA, nrow(best_train_pred_prob)),
                                Prediction_Prob =  train_prob,  # add prediction probability
                                OS_month = Train$OS_month,
                                OS_event = Train$OS)

    Test_results <- data.frame(Sample_ID = rownames(Test),  # Sample ID column
                               Actual = Test$Actual_Risk_Group,
                               Predicted_Risk_Group = best_test_pred,  # Add Predicted Risk Group
                               High_Risk_Prob = if ("High_Risk" %in% colnames(best_test_pred_prob)) best_test_pred_prob[, "High_Risk"] else rep(NA, nrow(best_test_pred_prob)),
                               Low_Risk_Prob = if ("Low_Risk" %in% colnames(best_test_pred_prob)) best_test_pred_prob[, "Low_Risk"] else rep(NA, nrow(best_test_pred_prob)),
                               Prediction_Prob = test_prob,  # add prediction probability
                               OS_month = Test$OS_month,
                               OS_event = Test$OS)


    row.names(Test_results) <- Test_results$Sample_ID
    row.names(Train_results) <- Train_results$Sample_ID


    # Extracting OOB Misclassification Rate
    OOB_misclassification <- best_model$err.rate[nrow(best_model$err.rate), 1]

    # Extracting class-wise misclassification errors
    Error_rates <- round(best_model$err.block.rate,3)
    All_class_error_rates <- Error_rates[1]
    #class_error_rates <- best_model$confusion[, "class.error"]
    high_risk_error <-  Error_rates[2]
    low_risk_error <-  Error_rates[3]

    # Compute confusion matrices
    cm_best_train  <- confusionMatrix(as.factor(best_train_pred), as.factor(Train$Actual_Risk_Group))
    cm_best_test <- confusionMatrix(as.factor(best_test_pred), as.factor(Test$Actual_Risk_Group))

    best_train_misclassification <- round(1 - cm_best_train$overall["Accuracy"], 3)
    best_test_misclassification <- round(1 - cm_best_test$overall["Accuracy"], 3)

    best_train_acc<- round(cm_best_train$overall["Accuracy"] * 100, 2)
    best_test_acc <- round(cm_best_test$overall["Accuracy"]* 100, 2)
    best_train_sens<- round(cm_best_train$byClass["Sensitivity"]* 100, 2)
    best_test_sens<- round(cm_best_test$byClass["Sensitivity"]* 100, 2)
    best_train_spec <- round(cm_best_train$byClass["Specificity"]* 100, 2)
    best_test_spec <- round(cm_best_test$byClass["Specificity"]* 100, 2)

    # Create a dataframe to store performance parameters for train and test data
    misclassification_results <- data.frame(
      Best_ntree = best_ntree,
      OOB_Misclassification = round(OOB_misclassification, 3),
      High_Risk_Error = round(high_risk_error, 3),
      Low_Risk_Error = round(low_risk_error, 3),
      Train_Misclassification_Error = best_train_misclassification,
      Train_Accuracy = best_train_acc,
      Train_Sensitivity = best_train_sens,
      Train_Specificity = best_train_spec,
      Test_Misclassification_Error = best_test_misclassification,
      Test_Accuracy = best_test_acc,
      Test_Sensitivity = best_test_sens,
      Test_Specificity = best_test_spec

    )


  }


  # Return results
  # Return results as a list
  return(list(
    Train_results = Train_results,
    Test_results = Test_results,
    best_model = best_model,
    All_tree_matrix =  all_metrics_matrix ,
    cm_best_train = cm_best_train,
    cm_best_test = cm_best_test,
    #best_test_pred = best_test_pred,
    misclassification_results = misclassification_results
  ))


}


