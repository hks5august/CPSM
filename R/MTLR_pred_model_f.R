#' @title MTLR Prediction Model Function
#'
#' @description This function builds a Multi-Task Logistic Regression (MTLR) survival model
#' based on clinical or selected feature or their combinations. It can
#' handle various model types based on user input (e.g., Clinical features
#' only, Clinical features with PI score, selected genes set).
#'
#' @name MTLR_pred_model_f
#' @param train_clin_data A data frame containing the clinical data for the
#' training set.
#' @param test_clin_data A data frame containing the clinical data for the
#' test set.
#' @param Model_type An integer indicating the model type:
#'   \itemize{
#'     \item 1: Clinical features model
#'     \item 2: PI score model
#'     \item 3: PI score with Clinical features model
#'      \item 4: selected gene set or clinical with selected gene set based model
#'   }
#' @param train_features_data A data frame containing the feature data for the
#' training set.
#' @param test_features_data A data frame containing the feature data for the
#' test set.
#' @param Clin_Feature_List A list of clinical features to be used in the
#' model.
#' @param surv_time A string specifying the name of the survival time column.
#' @param surv_event A string specifying the name of the survival event column.
#'
#' @return A list containing the following components:
#' \item{Model}{The fitted MTLR model.}
#' \item{Error_matrix}{A matrix of error evaluation metrics ( C-index and MAE)
#' for training and test sets.}
#' \item{Survival_predictions}{Predicted survival probabilities, mean survival
#' times, and median survival times.}
#' \item{Model_type}{The model type selected by the user.}
#'
#' @import reshape2
#' @import survival
#' @import survminer
#' @import MTLR
#' @importFrom stats median complete.cases
#'
#' @examples
#' # Example usage of the MTLR_pred_model_f function
#' data(Train_Clin, package = "CPSM")
#' data(Test_Clin, package = "CPSM")
#' data(Key_Clin_feature_list, package = "CPSM")
#' result <- MTLR_pred_model_f(
#'    train_clin_data = Train_Clin,
#'    test_clin_data = Test_Clin,
#'    Model_type = 1,
#'    train_features_data = Train_Clin,
#'    test_features_data = Test_Clin,
#'    Clin_Feature_List = Key_Clin_feature_list,
#'    surv_time = "OS_month",
#'    surv_event = "OS"
#'   )
#'
#'
#' @export



MTLR_pred_model_f <- function(train_clin_data, test_clin_data, Model_type,
                              train_features_data, test_features_data,
                              Clin_Feature_List, surv_time, surv_event) {
  
  # Validate input
  if (any(sapply(list(train_clin_data, test_clin_data, Model_type,
                      train_features_data, test_features_data,  
                      Clin_Feature_List, surv_time, surv_event), length) == 0)) {
    stop("Error: Empty input variable detected.")
  }
  
  # Prepare clinical data
  tr_clin <- train_clin_data
  te_clin <- test_clin_data
  colnames(tr_clin)[colnames(tr_clin) == surv_time] <- "OS_month"
  colnames(tr_clin)[colnames(tr_clin) == surv_event] <- "OS"
  colnames(te_clin)[colnames(te_clin) == surv_time] <- "OS_month"
  colnames(te_clin)[colnames(te_clin) == surv_event] <- "OS"
  
  # Combine clinical and features data
  train_all <- cbind(tr_clin, train_features_data)
  test_all <- cbind(te_clin, test_features_data)
  
  #rename  Clin_Feature_List to Feature_List
  Feature_List <- Clin_Feature_List
  # Feature list from Feature_List
  if ("ID" %in% colnames(Feature_List)) {
    # If "ID" is present as a column name
    features <- Feature_List$ID
  } else {
    # If "ID" is not found, check if there is another column with a different name
    # You can add more checks here if you expect the name to vary (e.g., "FeatureID", "GeneList", etc.)
    alternative_column <- colnames(Feature_List)[1]  # Assuming the first column might be the feature list
    message("The expected 'ID' column was not found. Using the first column: ", alternative_column)
    features <- Feature_List[[alternative_column]]
  }
  # Step 2: Identify missing features in train/test
  missing_features_tr <- setdiff(features, colnames(train_features_data))
  missing_features_te <- setdiff(features, colnames(test_features_data))
  
  # Step 3: Show warnings if features are missing
  if (length(missing_features_tr) > 0) {
    warning("The following features are missing in the training set: ", paste(missing_features_tr, collapse = ", "))
  }
  if (length(missing_features_te) > 0) {
    warning("The following features are missing in the test set: ", paste(missing_features_te, collapse = ", "))
  }
  
  # Step 4: Use only features common to both train and test
  common_features <- intersect(features, intersect(colnames(train_features_data), colnames(test_features_data)))
  
  if (length(common_features) == 0) {
    stop("Error: No common features found between training and test datasets.")
  } else {
    message("Using ", length(common_features), " common features for modeling: ", paste(common_features, collapse = ", "))
  }
  
  # Select model based on Model_type
  if (Model_type == 1) {
    # Clinical features only
    tr_data <- cbind(tr_clin[c("OS", "OS_month")], tr_clin[, common_features, drop = FALSE])
    te_data <- cbind(te_clin[c("OS", "OS_month")], te_clin[, common_features, drop = FALSE])
    
  } else if (Model_type == 2) {
    # PI only
    tr_data <- cbind(tr_clin[c("OS", "OS_month")], train_all["PI"])
    te_data <- cbind(te_clin[c("OS", "OS_month")], test_all["PI"])
    
  } else if (Model_type == 3) {
    # Clinical + PI
    tr_data <- cbind(tr_clin[c("OS", "OS_month")], train_all[, common_features, drop = FALSE])
    te_data <- cbind(te_clin[c("OS", "OS_month")], test_all[, common_features, drop = FALSE])
    
  } else if (Model_type == 4) {
    # Model with selected genes or combination of clinical features and genes
    tr_data <- cbind(tr_clin[c("OS", "OS_month")], train_all[, common_features, drop = FALSE])
    te_data <- cbind(te_clin[c("OS", "OS_month")], test_all[, common_features, drop = FALSE])
    
  } else {
    stop("Invalid Model_type. Must be 1, 2, 3, or 4.")
  }
  
  # Remove NA
  sel_clin_tr2 <- na.omit(tr_data)
  sel_clin_te2 <- na.omit(te_data)
  
  # Create MTLR model
  formula <- survival::Surv(OS_month, OS) ~ .
  # Make model
  Mod <- MTLR::mtlr(formula = formula, data = sel_clin_tr2)
  
  # Predictions on training data
  survival_curves_tr <- predict(Mod, sel_clin_tr2, type = "survivalcurve")
  mean_survival_tr <- predict(Mod, sel_clin_tr2, type = "mean_time")
  median_survival_tr <- predict(Mod, sel_clin_tr2, type = "median_time")
  event_probabilities_tr <- predict(Mod, sel_clin_tr2, type = "prob_event")
  
  # Prepare survival curve data
  colnames(survival_curves_tr) <- c("time_point", rownames(sel_clin_tr2))
  survCurves_tr_df <- as.data.frame(survival_curves_tr)
  
  survival_results_tr <- cbind(
    Mean_Survival = mean_survival_tr,
    Median_Survival = median_survival_tr,
    Event_Probability = event_probabilities_tr,
    Actual_OS_Time = sel_clin_tr2$OS_month,
    OS_Event = sel_clin_tr2$OS
  )
  rownames(survival_results_tr) <- rownames(sel_clin_tr2)
  
  # Predictions on test data
  survival_curves_te <- predict(Mod, sel_clin_te2, type = "survivalcurve")
  mean_survival_te <- predict(Mod, sel_clin_te2, type = "mean_time")
  median_survival_te <- predict(Mod, sel_clin_te2, type = "median_time")
  event_probabilities_te <- predict(Mod, sel_clin_te2, type = "prob_event")
  
  # Prepare survival curve data
  colnames(survival_curves_te) <- c("time_point", rownames(sel_clin_te2))
  survCurves_te_df <- as.data.frame(survival_curves_te)
  
  survival_results_te <- cbind(
    Mean_Survival = mean_survival_te,
    Median_Survival = median_survival_te,
    Event_Probability = event_probabilities_te,
    Actual_OS_Time = sel_clin_te2$OS_month,
    OS_Event = sel_clin_te2$OS
  )
  rownames(survival_results_te) <- rownames(sel_clin_te2)
  
  # Prepare mean/median survival summary
  survival_summary_te <- cbind(
    ID = rownames(sel_clin_te2),
    Mean = mean_survival_te,
    Median = median_survival_te,
    OS_month = sel_clin_te2$OS_month
  )
  
  # Prepare mean/median survival summary for training data
  survival_summary_tr <- cbind(
    ID = rownames(sel_clin_tr2),
    Mean = mean_survival_tr,
    Median = median_survival_tr,
    OS_month = sel_clin_tr2$OS_month
  )
  
  # C-Index and MAE calculation
  # C-Index calculation
  
  # Create survival object for training data
  surv_obj1_tr <- Surv(sel_clin_tr2$OS_month, sel_clin_tr2$OS)
  
  # Calculate C-index for training data
  c_index1_tr <- round(concordance(surv_obj1_tr ~ median_survival_tr)$concordance, 2)
  
  # Create survival object for test data
  surv_obj1_te <- Surv(sel_clin_te2$OS_month, sel_clin_te2$OS)
  
  # Calculate C-index for test data
  c_index1_te <- round(concordance(surv_obj1_te ~ median_survival_te)$concordance, 2)
  
  # --- MAE Calculation ----
  # Convert training survival summary to a data frame and ensure numeric values
  survival_summary_tr <- as.data.frame(survival_summary_tr)
  survival_summary_tr$OS_month <- as.numeric(survival_summary_tr$OS_month)
  survival_summary_tr$Mean <- as.numeric(survival_summary_tr$Mean)
  survival_summary_tr$Median <- as.numeric(survival_summary_tr$Median)
  
  # Compute MAE for training data
  mean_mae_tr <- round(mean(abs(survival_summary_tr$OS_month - survival_summary_tr$Mean), na.rm = TRUE), 2)
  median_mae_tr <- round(median(abs(survival_summary_tr$OS_month - survival_summary_tr$Median), na.rm = TRUE), 2)
  
  # Convert test survival summary to a data frame and ensure numeric values
  survival_summary_te <- as.data.frame(survival_summary_te)
  survival_summary_te$OS_month <- as.numeric(survival_summary_te$OS_month)
  survival_summary_te$Mean <- as.numeric(survival_summary_te$Mean)
  survival_summary_te$Median <- as.numeric(survival_summary_te$Median)
  
  # Compute MAE for test data
  mean_mae_te <- round(mean(abs(survival_summary_te$OS_month - survival_summary_te$Mean), na.rm = TRUE), 2)
  median_mae_te <- round(median(abs(survival_summary_te$OS_month - survival_summary_te$Median), na.rm = TRUE), 2)
  
  Error_mat_tr <- cbind(c_index1_tr,  mean_mae_tr,  median_mae_tr)
  Error_mat_te <- cbind(c_index1_te,  mean_mae_te,  median_mae_te)
  
  Error_mat <- rbind(Error_mat_tr , Error_mat_te)
  colnames(Error_mat) <- c("C_index", "Mean_MAE", "Median_MAE")
  rownames(Error_mat) <- c("Training_set", "Test_set")
  
  
  #prepare final outputs
  survCurves_df <- survCurves_te_df
  selected_train_data = as.data.frame(sel_clin_tr2)
  selected_test_data = as.data.frame(sel_clin_te2)
  mean_median_surv_d = survival_summary_te
  Error_mat = Error_mat
  surv_res = survival_results_te
  
  # Return results
  return(list(
    survCurves_data = survCurves_te_df,
    mean_median_survival_time_data = survival_summary_te,
    survival_result_based_on_MTLR = survival_results_te,
    Error_mat_for_Model = Error_mat,
    selected_train_data = selected_train_data,
    selected_test_data = selected_test_data
  ))
}



