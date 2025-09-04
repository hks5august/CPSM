#' MTLR Prediction Model Function
#'
#' This function builds a Multi-Task Logistic Regression (MTLR) survival model
#' based on clinical or selected feature or their combinations. It can
#' handle various model types based on user input (e.g., Clinical features
#' only, Clinical features with PI score, selected genes set).
#'
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
#' @importFrom SurvMetrics IBS
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
  # Check if any input variable is empty
  if (length(train_clin_data) == 0 || length(test_clin_data) == 0 ||
    length(Model_type) == 0 || length(train_features_data) == 0 ||
    length(test_features_data) == 0 || length(Clin_Feature_List) == 0 ||
    length(surv_time) == 0 || length(surv_event) == 0) {
    message("Error: Empty input variable detected.")
  }

  # load data
  tr_clin1 <- train_clin_data
  te_clin1 <- test_clin_data
  # rename survival time and event column name
  colnames(tr_clin1)[colnames(tr_clin1) == surv_time] <- "OS_month"
  colnames(tr_clin1)[colnames(tr_clin1) == surv_event] <- "OS"
  # test
  colnames(te_clin1)[colnames(te_clin1) == surv_time] <- "OS_month"
  colnames(te_clin1)[colnames(te_clin1) == surv_event] <- "OS"

  # load data
  train_features_data1 <- train_features_data
  test_features_data1 <- test_features_data
  # combine clinical and feature data
  tr_data2 <- cbind(tr_clin1, train_features_data1)
  te_data2 <- cbind(te_clin1, test_features_data1)
  # Load user defined a list of features forclinical data
  ftr_list <- Clin_Feature_List
  
# model1 - MTLR Model with Selected Clin features
  if (Model_type == 1) {
    # Load user defined a list of features for clin data
    ftr_list <- Clin_Feature_List
    # create data frame with selected features (user provided list)
    sel_clin_tr <- as.data.frame(tr_clin1[, colnames(tr_clin1) %in%
        c(ftr_list$ID), ])

    sel_clin_te <- as.data.frame(te_clin1[, colnames(te_clin1) %in%
        c(ftr_list$ID), ])
    # add survival information
    sel_clin_tr1 <- cbind(tr_clin1["OS"], tr_clin1["OS_month"], sel_clin_tr)
    sel_clin_te1 <- cbind(te_clin1["OS"], te_clin1["OS_month"], sel_clin_te)
    # create training and test data after removing NA values
    sel_clin_te2 <- na.omit(sel_clin_te1)
    sel_clin_tr2 <- na.omit(sel_clin_tr1)
    # create MTLR  model
    formula1 <- survival::Surv(OS_month, OS) ~ .
    # make  model!
    Mod1 <- MTLR::mtlr(formula = formula1, data = sel_clin_tr2)
    
    # Predictions on training data
    survival_curves_tr <- predict(Mod1, sel_clin_tr2, type = "survivalcurve")
    mean_survival_tr <- predict(Mod1, sel_clin_tr2, type = "mean_time")
    median_survival_tr <- predict(Mod1, sel_clin_tr2, type = "median_time")
    event_probabilities_tr <- predict(Mod1, sel_clin_tr2, type = "prob_event")
    
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
    survival_curves_te <- predict(Mod1, sel_clin_te2, type = "survivalcurve")
    mean_survival_te <- predict(Mod1, sel_clin_te2, type = "mean_time")
    median_survival_te <- predict(Mod1, sel_clin_te2, type = "median_time")
    event_probabilities_te <- predict(Mod1, sel_clin_te2, type = "prob_event")
    
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
    
    
    # Prepare mean/median survival summary
    survival_summary_tr <- cbind(
      ID = rownames(sel_clin_tr2),
      Mean = mean_survival_tr,
      Median = median_survival_tr,
      OS_month = sel_clin_tr2$OS_month
    )
    
    
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
   
    # IBS calculation
    # Training data
    # Survival probabilities at event times
    surv_probs_tr <- predict(Mod1, sel_clin_tr2, type = "prob_times")

    # Matrix of survival probabilities (drop time column)
    sp_matrix_tr <- as.matrix(surv_probs_tr[ , -1])

    #Extract the time grid for IBSrange
    time_points_tr <- surv_probs_tr$time[-1]

    #Integrated Brier Score with survmetrics
    ibs_tr <- IBS(
    object   = surv_obj1_tr,
    sp_matrix = sp_matrix_tr,
    IBSrange  = time_points_tr)

   # Round up value
    ibs_tr <- round(ibs_tr, 3)

    #IBS calculation for Test data
    # Predicted survival probabilities at event times
    surv_probs_te <- predict(Mod1, sel_clin_te2, type = "prob_times")

    # Matrix of survival probabilities (drop time column)
    sp_matrix_te <- as.matrix(surv_probs_te[ , -1])

    #Extract the time grid for IBSrange
    time_points_te <- surv_probs_te$time[-1]

    #Integrated Brier Score with survmetrics
    ibs_te <- IBS(
    object   = surv_obj1_te,
    sp_matrix = sp_matrix_te,
    IBSrange  = time_points_te)

  # Round up value
  ibs_te <- round(ibs_te, 3) 

    Error_mat_tr <- cbind(c_index1_tr, mean_mae_tr, median_mae_tr, round(ibs_tr, 3))
    Error_mat_te <- cbind(c_index1_te, mean_mae_te, median_mae_te, round(ibs_te, 3))

    colnames(Error_mat_tr) <- colnames(Error_mat_te) <- c("C_index", "Mean_MAE", "Median_MAE", "IBS")
    rownames(Error_mat) <- c("Training_set", "Test_set")

 }

# Model2 -  Model with only PI score

  else if (Model_type == 2) {
    # combine clinical and feature data
    sel_clin_tr1 <- cbind(
      tr_clin1["OS"], tr_clin1["OS_month"],
      tr_data2["PI"]
    )
    sel_clin_te1 <- cbind(
      te_clin1["OS"], te_clin1["OS_month"],
      te_data2["PI"]
    )
    # create training and test data after removing NA values
    sel_clin_tr2 <- na.omit(sel_clin_tr1)
    sel_clin_te2 <- na.omit(sel_clin_te1)

    # create MTLR  model
    formula2 <- survival::Surv(OS_month, OS) ~ .

    Mod2 <- mtlr(formula = formula2, data = sel_clin_tr2)

    # Model Predictions
    # Predictions on training data
    survival_curves_tr2 <- predict(Mod2, sel_clin_tr2, type = "survivalcurve")
    mean_survival_tr2 <- predict(Mod2, sel_clin_tr2, type = "mean_time")
    median_survival_tr2 <- predict(Mod2, sel_clin_tr2, type = "median_time")
    event_probabilities_tr2 <- predict(Mod2, sel_clin_tr2, type = "prob_event")
    
    # Prepare survival curve data       
    colnames(survival_curves_tr2) <- c("time_point", rownames(sel_clin_tr2))
    survCurves_tr_df2 <- as.data.frame(survival_curves_tr2)
    
    
    survival_results_tr2 <- cbind(
      Mean_Survival = mean_survival_tr2,
      Median_Survival = median_survival_tr2,
      Event_Probability = event_probabilities_tr2,
      Actual_OS_Time = sel_clin_tr2$OS_month,
      OS_Event = sel_clin_tr2$OS
    )
    rownames(survival_results_tr2) <- rownames(sel_clin_tr2)
    
    # Predictions on test data
    survival_curves_te2 <- predict(Mod2, sel_clin_te2, type = "survivalcurve")
    mean_survival_te2 <- predict(Mod2, sel_clin_te2, type = "mean_time")
    median_survival_te2 <- predict(Mod2, sel_clin_te2, type = "median_time")
    event_probabilities_te2 <- predict(Mod2, sel_clin_te2, type = "prob_event")
    
    # Prepare survival curve data       
    colnames(survival_curves_te2) <- c("time_point", rownames(sel_clin_te2))
    survCurves_te_df2 <- as.data.frame(survival_curves_te2)
    
    survival_results_te2 <- cbind(
      Mean_Survival = mean_survival_te2,
      Median_Survival = median_survival_te2,
      Event_Probability = event_probabilities_te2,
      Actual_OS_Time = sel_clin_te2$OS_month,
      OS_Event = sel_clin_te2$OS
    )
    rownames(survival_results_te2) <- rownames(sel_clin_te2)
    
    
    # Prepare mean/median survival summary
    survival_summary_te2 <- cbind(
      ID = rownames(sel_clin_te2),
      Mean = mean_survival_te2,
      Median = median_survival_te2,
      OS_month = sel_clin_te2$OS_month
    )
    
    
    # Prepare mean/median survival summary
    survival_summary_tr2 <- cbind(
      ID = rownames(sel_clin_tr2),
      Mean = mean_survival_tr2,
      Median = median_survival_tr2,
      OS_month = sel_clin_tr2$OS_month
    )
    
    
    # C-Index calculation
    
    # Create survival object for training data
    surv_obj1_tr2 <- Surv(sel_clin_tr2$OS_month, sel_clin_tr2$OS)
    
    # Calculate C-index for training data
    c_index1_tr2 <- round(concordance(surv_obj1_tr2 ~ median_survival_tr2)$concordance, 2)
    
    # Create survival object for test data
    surv_obj1_te2 <- Surv(sel_clin_te2$OS_month, sel_clin_te2$OS)
    
    # Calculate C-index for test data
    c_index1_te2 <- round(concordance(surv_obj1_te2 ~ median_survival_te2)$concordance, 2)
    
    # --- MAE Calculation ----
    # Convert training survival summary to a data frame and ensure numeric values
    survival_summary_tr2 <- as.data.frame(survival_summary_tr2)
    survival_summary_tr2$OS_month <- as.numeric(survival_summary_tr2$OS_month)
    survival_summary_tr2$Mean <- as.numeric(survival_summary_tr2$Mean)
    survival_summary_tr2$Median <- as.numeric(survival_summary_tr2$Median)
    
    # Compute MAE for training data
    mean_mae_tr2 <- round(mean(abs(survival_summary_tr2$OS_month - survival_summary_tr2$Mean), na.rm = TRUE), 2)
    median_mae_tr2 <- round(median(abs(survival_summary_tr2$OS_month - survival_summary_tr2$Median), na.rm = TRUE), 2)
    
    # Convert test survival summary to a data frame and ensure numeric values
    survival_summary_te2 <- as.data.frame(survival_summary_te2)
    survival_summary_te2$OS_month <- as.numeric(survival_summary_te2$OS_month)
    survival_summary_te2$Mean <- as.numeric(survival_summary_te2$Mean)
    survival_summary_te2$Median <- as.numeric(survival_summary_te2$Median)
    
    # Compute MAE for test data
    mean_mae_te2 <- round(mean(abs(survival_summary_te2$OS_month - survival_summary_te2$Mean), na.rm = TRUE), 2)
    median_mae_te2 <- round(median(abs(survival_summary_te2$OS_month - survival_summary_te2$Median), na.rm = TRUE), 2)
    
    Error_mat_tr2 <- cbind(c_index1_tr2,  mean_mae_tr2,  median_mae_tr2)
    Error_mat_te2 <- cbind(c_index1_te2,  mean_mae_te2,  median_mae_te2)
    
    Error_mat2 <- rbind(Error_mat_tr2 , Error_mat_te2)
    colnames(Error_mat2) <- c("C_index", "Mean_MAE", "Median_MAE")
    rownames(Error_mat2) <- c("Training_set", "Test_set")
    
    #IBS calculation
    
#model3  
} else if (Model_type == 3) { # Model3- Model with PI & Clin features
    # create data frame with selected features (user provided list)
    sel_clin_tr <- as.data.frame(tr_data2[, colnames(tr_data2) %in%
        c(ftr_list$ID), ])
    sel_clin_te <- as.data.frame(te_data2[, colnames(te_data2) %in%
        c(ftr_list$ID), ])

    # add survival information
    sel_clin_tr1 <- cbind(tr_clin1["OS"], tr_clin1["OS_month"], sel_clin_tr)
    sel_clin_te1 <- cbind(te_clin1["OS"], te_clin1["OS_month"], sel_clin_te)

    # remove samples where info missing for any of selected feature in
    # training data
    sel_clin_tr2 <- na.omit(sel_clin_tr1)
    # test data
    sel_clin_te2 <- na.omit(sel_clin_te1)
    # develop MTLR model
    # create formula
    formula3 <- survival::Surv(OS_month, OS) ~ .
    #  make our model
    Mod3 <- mtlr(formula = formula3, data = sel_clin_tr2)

    # Model Predictions
    # Predictions on training data
    survival_curves_tr3 <- predict(Mod3, sel_clin_tr2, type = "survivalcurve")
    mean_survival_tr3 <- predict(Mod3, sel_clin_tr2, type = "mean_time")
    median_survival_tr3 <- predict(Mod3, sel_clin_tr2, type = "median_time")
    event_probabilities_tr3 <- predict(Mod3, sel_clin_tr2, type = "prob_event")
    
    # Prepare survival curve data       
    colnames(survival_curves_tr3) <- c("time_point", rownames(sel_clin_tr2))
    survCurves_tr_df3 <- as.data.frame(survival_curves_tr3)
    
    
    survival_results_tr3 <- cbind(
      Mean_Survival = mean_survival_tr3,
      Median_Survival = median_survival_tr3,
      Event_Probability = event_probabilities_tr3,
      Actual_OS_Time = sel_clin_tr2$OS_month,
      OS_Event = sel_clin_tr2$OS
    )
    rownames(survival_results_tr3) <- rownames(sel_clin_tr2)
    
    # Predictions on test data
    survival_curves_te3 <- predict(Mod3, sel_clin_te2, type = "survivalcurve")
    mean_survival_te3 <- predict(Mod3, sel_clin_te2, type = "mean_time")
    median_survival_te3 <- predict(Mod3, sel_clin_te2, type = "median_time")
    event_probabilities_te3 <- predict(Mod3, sel_clin_te2, type = "prob_event")
    
    # Prepare survival curve data       
    colnames(survival_curves_te3) <- c("time_point", rownames(sel_clin_te2))
    survCurves_te_df3 <- as.data.frame(survival_curves_te3)
    
    survival_results_te3 <- cbind(
      Mean_Survival = mean_survival_te3,
      Median_Survival = median_survival_te3,
      Event_Probability = event_probabilities_te3,
      Actual_OS_Time = sel_clin_te2$OS_month,
      OS_Event = sel_clin_te2$OS
    )
    rownames(survival_results_te3) <- rownames(sel_clin_te2)
    
    
    # Prepare mean/median survival summary
    survival_summary_te3 <- cbind(
      ID = rownames(sel_clin_te2),
      Mean = mean_survival_te3,
      Median = median_survival_te3,
      OS_month = sel_clin_te2$OS_month
    )
    
    
    # Prepare mean/median survival summary
    survival_summary_tr3 <- cbind(
      ID = rownames(sel_clin_tr2),
      Mean = mean_survival_tr3,
      Median = median_survival_tr3,
      OS_month = sel_clin_tr2$OS_month
    )
    
    
    # C-Index calculation
    
    # Create survival object for training data
    surv_obj1_tr3 <- Surv(sel_clin_tr2$OS_month, sel_clin_tr2$OS)
    
    # Calculate C-index for training data
    c_index1_tr3 <- round(concordance(surv_obj1_tr3 ~ median_survival_tr3)$concordance, 2)
    
    # Create survival object for test data
    surv_obj1_te3 <- Surv(sel_clin_te2$OS_month, sel_clin_te2$OS)
    
    # Calculate C-index for test data
    c_index1_te3 <- round(concordance(surv_obj1_te3 ~ median_survival_te3)$concordance, 2)
    
    # --- MAE Calculation ----
    # Convert training survival summary to a data frame and ensure numeric values
    survival_summary_tr3 <- as.data.frame(survival_summary_tr3)
    survival_summary_tr3$OS_month <- as.numeric(survival_summary_tr3$OS_month)
    survival_summary_tr3$Mean <- as.numeric(survival_summary_tr3$Mean)
    survival_summary_tr3Median <- as.numeric(survival_summary_tr3$Median)
    
    # Compute MAE for training data
    mean_mae_tr3 <- round(mean(abs(survival_summary_tr3$OS_month - survival_summary_tr3$Mean), na.rm = TRUE), 2)
    median_mae_tr3 <- round(median(abs(survival_summary_tr3$OS_month - survival_summary_tr3$Median), na.rm = TRUE), 2)
    
    # Convert test survival summary to a data frame and ensure numeric values
    survival_summary_te3 <- as.data.frame(survival_summary_te3)
    survival_summary_te3$OS_month <- as.numeric(survival_summary_te3$OS_month)
    survival_summary_te3$Mean <- as.numeric(survival_summary_te3$Mean)
    survival_summary_te3$Median <- as.numeric(survival_summary_te3$Median)
    
    # Compute MAE for test data
    mean_mae_te3 <- round(mean(abs(survival_summary_te3$OS_month - survival_summary_te3$Mean), na.rm = TRUE), 2)
    median_mae_te3 <- round(median(abs(survival_summary_te3$OS_month - survival_summary_te3$Median), na.rm = TRUE), 2)
    
    Error_mat_tr3 <- cbind(c_index1_tr3,  mean_mae_tr3,  median_mae_tr3)
    Error_mat_te3 <- cbind(c_index1_te3,  mean_mae_te3,  median_mae_te3)
    
    Error_mat3 <- rbind(Error_mat_tr3 , Error_mat_te3)
    colnames(Error_mat3) <- c("C_index", "Mean_MAE", "Median_MAE")
    rownames(Error_mat3) <- c("Training_set", "Test_set")
    
    # IBS calculation
    
    
  } else if (Model_type == 4) {
    ## Univariate with Clin features ##
    # create data frame with selected features (user provided list)
    sel_clin_tr <- as.data.frame(tr_data2[, colnames(tr_data2) %in%
      c(ftr_list$ID), ])
    sel_clin_te <- as.data.frame(te_data2[, colnames(te_data2) %in%
      c(ftr_list$ID), ])

    # add survival information
    sel_clin_tr1 <- cbind(tr_clin1["OS"], tr_clin1["OS_month"], sel_clin_tr)
    sel_clin_te1 <- cbind(te_clin1["OS"], te_clin1["OS_month"], sel_clin_te)

    # remove samples where features used in the models are missing
    sel_clin_tr2 <- na.omit(sel_clin_tr1)
    sel_clin_te2 <- na.omit(sel_clin_te1)
    # develop model using MTLR
    # craete formula
    formula5 <- survival::Surv(OS_month, OS) ~ .
    # Next, we just need the data argument which in our case is training.
    # We can finally make our first model!
    Mod5 <- mtlr(formula = formula5, data = sel_clin_tr2)

    # Model Predictions
    # Predictions on training data
    survival_curves_tr5 <- predict(Mod5, sel_clin_tr2, type = "survivalcurve")
    mean_survival_tr5 <- predict(Mod5, sel_clin_tr2, type = "mean_time")
    median_survival_tr5 <- predict(Mod5, sel_clin_tr2, type = "median_time")
    event_probabilities_tr5 <- predict(Mod5, sel_clin_tr2, type = "prob_event")
    
    # Prepare survival curve data       
    colnames(survival_curves_tr5) <- c("time_point", rownames(sel_clin_tr2))
    survCurves_tr_df5 <- as.data.frame(survival_curves_tr5)
    
    
    survival_results_tr5 <- cbind(
      Mean_Survival = mean_survival_tr5,
      Median_Survival = median_survival_tr5,
      Event_Probability = event_probabilities_tr5,
      Actual_OS_Time = sel_clin_tr2$OS_month,
      OS_Event = sel_clin_tr2$OS
    )
    rownames(survival_results_tr5) <- rownames(sel_clin_tr2)
    
    # Predictions on test data
    survival_curves_te5 <- predict(Mod5, sel_clin_te2, type = "survivalcurve")
    mean_survival_te5 <- predict(Mod5, sel_clin_te2, type = "mean_time")
    median_survival_te5 <- predict(Mod5, sel_clin_te2, type = "median_time")
    event_probabilities_te5 <- predict(Mod5, sel_clin_te2, type = "prob_event")
    
    # Prepare survival curve data       
    colnames(survival_curves_te5) <- c("time_point", rownames(sel_clin_te2))
    survCurves_te_df5 <- as.data.frame(survival_curves_te5)
    
    survival_results_te5 <- cbind(
      Mean_Survival = mean_survival_te5,
      Median_Survival = median_survival_te5,
      Event_Probability = event_probabilities_te5,
      Actual_OS_Time = sel_clin_te2$OS_month,
      OS_Event = sel_clin_te2$OS
    )
    rownames(survival_results_te5) <- rownames(sel_clin_te2)
    
    
    # Prepare mean/median survival summary
    survival_summary_te5 <- cbind(
      ID = rownames(sel_clin_te2),
      Mean = mean_survival_te5,
      Median = median_survival_te5,
      OS_month = sel_clin_te2$OS_month
    )
    
    
    # Prepare mean/median survival summary
    survival_summary_tr5 <- cbind(
      ID = rownames(sel_clin_tr2),
      Mean = mean_survival_tr5,
      Median = median_survival_tr5,
      OS_month = sel_clin_tr2$OS_month
    )
    
    
    # C-Index calculation
    
    # Create survival object for training data
    surv_obj1_tr5 <- Surv(sel_clin_tr2$OS_month, sel_clin_tr2$OS)
    
    # Calculate C-index for training data
    c_index1_tr5 <- round(concordance(surv_obj1_tr5 ~ median_survival_tr5)$concordance, 2)
    
    # Create survival object for test data
    surv_obj1_te5 <- Surv(sel_clin_te2$OS_month, sel_clin_te2$OS)
    
    # Calculate C-index for test data
    c_index1_te5 <- round(concordance(surv_obj1_te5 ~ median_survival_te5)$concordance, 2)
    
    # --- MAE Calculation ----
    # Convert training survival summary to a data frame and ensure numeric values
    survival_summary_tr5 <- as.data.frame(survival_summary_tr5)
    survival_summary_tr5$OS_month <- as.numeric(survival_summary_tr5$OS_month)
    survival_summary_tr5$Mean <- as.numeric(survival_summary_tr5$Mean)
    survival_summary_tr5Median <- as.numeric(survival_summary_tr5$Median)
    
    # Compute MAE for training data
    mean_mae_tr5 <- round(mean(abs(survival_summary_tr5$OS_month - survival_summary_tr5$Mean), na.rm = TRUE), 2)
    median_mae_tr5 <- round(mean(abs(survival_summary_tr5$OS_month -  as.numeric(survival_summary_tr5$Median)), na.rm = TRUE), 2)
    
    # Convert test survival summary to a data frame and ensure numeric values
    survival_summary_te5 <- as.data.frame(survival_summary_te5)
    survival_summary_te5$OS_month <- as.numeric(survival_summary_te5$OS_month)
    survival_summary_te5$Mean <- as.numeric(survival_summary_te5$Mean)
    survival_summary_te5$Median <- as.numeric(survival_summary_te5$Median)
    
    # Compute MAE for test data
    mean_mae_te5 <- round(mean(abs(survival_summary_te5$OS_month - survival_summary_te5$Mean), na.rm = TRUE), 2)
    median_mae_te5 <- round(median(abs(survival_summary_te5$OS_month - survival_summary_te5$Median), na.rm = TRUE), 2)
    
    Error_mat_tr5 <- cbind(c_index1_tr5,  mean_mae_tr5,  median_mae_tr5)
    Error_mat_te5 <- cbind(c_index1_te5,  mean_mae_te5,  median_mae_te5)
    
    Error_mat5 <- rbind(Error_mat_tr5 , Error_mat_te5)
    colnames(Error_mat5) <- c("C_index", "Mean_MAE", "Median_MAE")
    rownames(Error_mat5) <- c("Training_set", "Test_set")
   
    #IBS calculation 
    
  }
  if (Model_type == 1) {
    survCurves_df <- survCurves_te_df 
    selected_train_data = as.data.frame(sel_clin_tr2)
    selected_test_data = as.data.frame(sel_clin_te2)
    mean_median_surv_d <- survival_summary_te
    Error_mat <- Error_mat
    surv_res <- survival_results_te
  }
  if (Model_type == 2) {
    survCurves_df <- survCurves_te_df2
    selected_train_data = as.data.frame(sel_clin_tr2)
    selected_test_data = as.data.frame(sel_clin_te2)
    mean_median_surv_d <- survival_summary_te2
    Error_mat <- Error_mat2
    surv_res <- survival_results_te2
  }
  if (Model_type == 3) {
    survCurves_df <- survCurves_te_df3
    selected_train_data = as.data.frame(sel_clin_tr2)
    selected_test_data = as.data.frame(sel_clin_te2)
    mean_median_surv_d <- survival_summary_te3
    Error_mat <- Error_mat3
    surv_res <- survival_results_te3
  }
  if (Model_type == 4) {
    survCurves_df <- survCurves_te_df5
    selected_train_data = as.data.frame(sel_clin_tr2)
    selected_test_data = as.data.frame(sel_clin_te2)
    mean_median_surv_d <- survival_summary_te5
    Error_mat <- Error_mat5
    surv_res <- survival_results_te5
  }
  
  # Return a list containing data.
  return(list(
    survCurves_data = survCurves_df,
    mean_median_survival_time_data = mean_median_surv_d,
    survival_result_based_on_MTLR = surv_res,
    Error_mat_for_Model = Error_mat,
    selected_train_data = selected_train_data, 
    selected_test_data = selected_test_data
    
  ))
}


