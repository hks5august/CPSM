#' This function can generate 5 Prediction models types based on :(1)Clinical
#' features,  (2) PI score , (3) PI score + Clin, (4) Significant Univariate
#' features, (5) Significant Univariate features + clin features, using MTLR
#' method.
#' Further, this model will predict the survival of test data in terms of
#' survival probability over different time points, mean and median survival
#' time pf patients in test data.

#' @param train_clin_data :args1 - training data with Clin features (Patients
#' data with clinical and gene expression, where samples are in rows and
#' features/genes are in columns)
#' @param test_clin_data :args2 -  test data with Clin features  (Patients data
#' with clinical and gene expression, where samples are in rows and
#' features/genes are in columns)
#' @param  Model_type: args3 - Prediction Model type : 1 - Clinical features,
#' 2 - PI score features, 3-  PI score + Clinical features ; 4 - Univariate
#' features, 5 -  Univariate features + Clinical features
#' @param train_features_data :args4 - Training data with selected features
#' only (Significant Univariate/LASSO PI score)
#' @param test_features_data :args5 - Training data with selected features
#' only (Significant Univariate/LASSO PI score)
#' @param Clin_Feature_List :args6 - List of key clinical features which
#' contain non-unique values (e.g. Key_Clin_feature_list.txt, Key_PI_list.txt
#' Key_Clin_features_with_PI_list.txt, Key_univariate_features_list.txt,
#' Key_univariate_features_with_Clin_list.txt)
#' @param surv_time :arg7 - name of column which contain survival time
#' (in days) information
#' @param surv_event :arg8 - name of column which contain survival
#' eventinformation
#' @import dplyr
#' @import reshape2
#' @import survival
#' @import survminer
#' @import MTLR
#' @import SurvMetrics
#' @import pec
#' @import ggplot2
#' @examples
#' data(Train_Clin, package = "CPSM")
#' data(Test_Clin, package = "CPSM")
#' data(Key_Clin_feature_list, package = "CPSM")
#' MTLR_pred_model_f(
#'   train_clin_data = Train_Clin, test_clin_data =
#'     test_clin_data, Model_type = 1, train_features_data = Train_Clin,
#'   test_features_data = test_clin_data, Clin_Feature_List =
#'     Key_Clin_feature_list, surv_time = "OS_month", surv_event = "OS"
#' )
#' Usage:MTLR_pred_model_f(
#'   train_clin_data, test_clin_data, Model_type,
#'   train_features_data, test_features_data, Clin_Feature_List, surv_time,
#'   surv_event
#' )
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
    # Prediction on Test Data
    survCurves1 <- predict(Mod1, sel_clin_te1, type = "survivalcurve")
    # define column names
    colnames(survCurves1) <- c("time_point", rownames(sel_clin_te2))
    survCurves1_df <- as.data.frame(survCurves1)
    survivalcurve1 <- predict(Mod1, sel_clin_te2, type = "survivalcurve")
    Survival_prob_event1 <- predict(Mod1, sel_clin_te2, type = "prob_event")
    # Predicted Mean
    meanSurv1_tr <- predict(Mod1, sel_clin_tr2, type = "mean_time")
    meanSurv1 <- predict(Mod1, sel_clin_te2, type = "mean_time")
    # Predicted Median
    medianSurv1_tr <- predict(Mod1, sel_clin_tr2, type = "median_time")
    medianSurv1 <- predict(Mod1, sel_clin_te2, type = "median_time")
    # create a dataframe of  predicted  mean survival time
    meanSurv_d1 <- as.data.frame(meanSurv1)
    # create a dataframe of  predicted  median survival time
    medianSurv_d1 <- as.data.frame(medianSurv1)
    names1 <- as.data.frame(rownames(sel_clin_te2))
    # create dataframe combining both predicted mean & median survivaltime
    mean_median_surv1_d <- cbind(names1, meanSurv_d1, medianSurv_d1)
    rownames(mean_median_surv1_d) <- c(rownames(sel_clin_te2))
    colnames(mean_median_surv1_d) <- c("IDs", "Mean", "Median")
    # Survival Probability at Event Time
    survivalProbs_p1_tr <- predict(Mod1, sel_clin_tr2, type = "prob_times")
    # extract prob times at diff times points
    survivalProbs_t_mat_tr <- as.matrix(survivalProbs_p1_tr)
    survivalProbs_t_mat1_tr <- survivalProbs_t_mat_tr[, -1]
    survivalProbs_t_mat1_t_tr <- t(survivalProbs_t_mat1_tr)
    survivalProbs_t_mat1_t2_tr <- survivalProbs_t_mat1_t_tr[, -1]
    # Test data
    survivalProbs_p1 <- predict(Mod1, sel_clin_te2, type = "prob_times")
    # extract prob times at diff times points
    survivalProbs_t_mat <- as.matrix(survivalProbs_p1)
    survivalProbs_t_mat1 <- survivalProbs_t_mat[, -1]
    survivalProbs_t_mat1_t <- t(survivalProbs_t_mat1)
    survivalProbs_t_mat1_t2 <- survivalProbs_t_mat1_t[, -1]
    # combine predicted mean, median, survival
    # probability and actual time and event
    surv_res1 <- cbind(
      meanSurv1, medianSurv1, Survival_prob_event1,
      sel_clin_te2$OS_month, sel_clin_te2$OS
    )
    # add column names and row names
    colnames(surv_res1) <- c(
      "Mean_Surv", "Median_surv", "Prob_event",
      "Actual_OS_time", "OS_event"
    )
    rownames(surv_res1) <- rownames(sel_clin_te2)
    # Calcualte Evalulation parameters on training data
    # create survival object
    surv_obj1_tr <- survival::Surv(sel_clin_tr2$OS_month, sel_clin_tr2$OS)
    # calculate IBS (Integrated Brier Score for test data
    IBS_1_tr <- round(IBS(surv_obj1_tr,
                          sp_matrix = survivalProbs_t_mat1_t2_tr,
                          survivalProbs_p1_tr$time[-1]), 3)

    # Calculate Concordance Index
    c_index1_tr <- round(SurvMetrics::Cindex(surv_obj1_tr,
                                             predicted = medianSurv1_tr), 2)
    # Combine evaluation parameters to get Matrix
    Error_mat_1_tr <- cbind(IBS_1_tr, c_index1_tr)
    # Calcualte Evalulation parameters on test data
    # create survival object
    surv_obj1 <- survival::Surv(sel_clin_te2$OS_month, sel_clin_te2$OS)

    # calculate IBS (Integrated Brier Score for test data
    IBS_1 <- round(IBS(surv_obj1,
                       sp_matrix = survivalProbs_t_mat1_t2,
                       survivalProbs_p1$time[-1]), 3)

    # Calculate Concordance Index
    c_index1 <- round(SurvMetrics::Cindex(surv_obj1,
                                          predicted = medianSurv1), 2)
    # Combine evaluation parameters to get Matrix
    Error_mat_1_te <- cbind(IBS_1, c_index1)
    Error_mat_1 <- rbind(Error_mat_1_tr, Error_mat_1_te)
    colnames(Error_mat_1) <- c("IBS_score", "c_index")
    rownames(Error_mat_1) <- c("Training_set", "Test_set")
  } else if (Model_type == 2) { # Model2- Model with only PI score
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
    # prediction on Test data
    survCurves2 <- predict(Mod2, sel_clin_te2, type = "survivalcurve")
    # add column names
    colnames(survCurves2) <- c("time_point", rownames(sel_clin_te2))

    # create dataframe of survival curve data
    survCurves2_df <- as.data.frame(survCurves2)

    # survivalcurve
    survivalcurve2 <- predict(Mod2, sel_clin_te2, type = "survivalcurve")
    # Prob_Event
    Survival_prob_event2_tr <- predict(Mod2, sel_clin_tr2,
      type = "prob_event"
    )
    Survival_prob_event2 <- predict(Mod2, sel_clin_te2, type = "prob_event")

    # Mean/Median Survival Time: In addition to the entire survival curve one
    # may also be interested in the average survival time.
    # This is again available from the predict function.
    # Predicted Mean
    meanSurv2_tr <- predict(Mod2, sel_clin_tr2, type = "mean_time")
    meanSurv2 <- predict(Mod2, sel_clin_te2, type = "mean_time")
    # Predicted Median
    medianSurv2_tr <- predict(Mod2, sel_clin_tr2, type = "median_time")
    medianSurv2 <- predict(Mod2, sel_clin_te2, type = "median_time")
    ## create a dataframe of  predicted  mean and median survival time
    meanSurv2_d <- as.data.frame(meanSurv2)
    medianSurv2_d <- as.data.frame(medianSurv2)

    # add rownames
    names_2 <- as.data.frame(rownames(sel_clin_te2))

    # create a dataframe combining both predicted  mean and median survival
    # time
    mean_median_surv2_d <- cbind(names_2, meanSurv2_d, medianSurv2_d)
    rownames(mean_median_surv2_d) <- c(rownames(sel_clin_te2))

    # add column names
    colnames(mean_median_surv2_d) <- c("IDs", "Mean", "Median")

    # Survival Probability at Event Time
    # training data prediction
    survivalProbs_p2_tr <- predict(Mod2, sel_clin_tr2, type = "prob_times")
    # extract prob times at diff times points
    survivalProbs_t_mat_2_tr <- as.matrix(survivalProbs_p2_tr)
    survivalProbs_t_mat1_2_tr <- survivalProbs_t_mat_2_tr[, -1]
    survivalProbs_t_mat1_t_2_tr <- t(survivalProbs_t_mat1_2_tr)
    survivalProbs_t_mat1_t2_2_tr <- survivalProbs_t_mat1_t_2_tr[, -1]
    # test data prediction
    survivalProbs_p2 <- predict(Mod2, sel_clin_te2, type = "prob_times")
    # extract prob times at diff times points
    survivalProbs_t_mat_2 <- as.matrix(survivalProbs_p2)
    survivalProbs_t_mat1_2 <- survivalProbs_t_mat_2[, -1]
    survivalProbs_t_mat1_t_2 <- t(survivalProbs_t_mat1_2)
    survivalProbs_t_mat1_t2_2 <- survivalProbs_t_mat1_t_2[, -1]

    # create a data frame combining predicted mean, median, survival
    # probability and actual time and event
    surv_res2 <- cbind(
      meanSurv2, medianSurv2, Survival_prob_event2,
      sel_clin_te2$OS_month, sel_clin_te2$OS
    )
    # add column and rownames
    colnames(surv_res2) <- c(
      "Mean_Surv", "Prob_event", "Median_surv",
      "Actual_OS_time", "OS_event"
    )
    rownames(surv_res2) <- rownames(sel_clin_te2)

    # Calcualte Evalulation/prediction parameters on training data
    # create survival object
    surv_obj_2_tr <- survival::Surv(sel_clin_tr2$OS_month, sel_clin_tr2$OS)

    # calculate IBS (Integrated Brier Score for test data
    IBS1_2_tr <- round(IBS(surv_obj_2_tr,
                           sp_matrix = survivalProbs_t_mat1_t2_2_tr,
                           survivalProbs_p2_tr$time[-1]), 3)

    # Concordance Index
    c_index_2_tr <- round(SurvMetrics::Cindex(surv_obj_2_tr,
                                              predicted = medianSurv2_tr), 2)

    # Combine evaluation parameters to get Matrix
    Error_mat_2_tr <- cbind(IBS1_2_tr, c_index_2_tr)

    # Calcualte Evalulation/prediction parameters on test data
    # create survival object
    surv_obj_2 <- survival::Surv(sel_clin_te2$OS_month, sel_clin_te2$OS)

    # calculate IBS (Integrated Brier Score for test data
    IBS1_2 <- round(IBS(surv_obj_2,
                        sp_matrix = survivalProbs_t_mat1_t2_2,
                        survivalProbs_p2$time[-1]), 3)

    # Concordance Index
    c_index_2 <- round(SurvMetrics::Cindex(surv_obj_2,
                                           predicted = medianSurv2), 2)

    # Combine evaluation parameters to get Matrix
    Error_mat_2_te <- cbind(IBS1_2, c_index_2)

    Error_mat_2 <- rbind(Error_mat_2_tr, Error_mat_2_te)
    colnames(Error_mat_2) <- c("IBS_Score", "c_index")
    rownames(Error_mat_2) <- c("Training_set", "Test_set")
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

    # survival curve data
    survCurves3 <- predict(Mod3, sel_clin_te2, type = "survivalcurve")
    # add column names
    colnames(survCurves3) <- c("time_point", rownames(sel_clin_te2))

    # make dataframe
    survCurves3_df <- as.data.frame(survCurves3)
    # Predicted Mean Survival Time
    # training data
    meanSurv3_tr <- predict(Mod3, sel_clin_tr2, type = "mean_time")
    # test data
    meanSurv3 <- predict(Mod3, sel_clin_te2, type = "mean_time")

    # Predicted Median survival time
    # training
    medianSurv3_tr <- predict(Mod3, sel_clin_tr2, type = "median_time")
    # test
    medianSurv3 <- predict(Mod3, sel_clin_te2, type = "median_time")
    # create dataframes of predicted mean and median sirvival time
    meanSurv3_d <- as.data.frame(meanSurv3)
    medianSurv3_d <- as.data.frame(medianSurv3)
    # combine both mean and median sirvival time
    # add rownames
    names_3 <- as.data.frame(rownames(sel_clin_te2))
    # combine both mean and median sirvival time
    mean_median_surv3_d <- cbind(names_3, meanSurv3_d, medianSurv3_d)

    # add row and column names
    rownames(mean_median_surv3_d) <- c(rownames(sel_clin_te2))
    colnames(mean_median_surv3_d) <- c("IDs", "Mean", "Median")

    # Survival Probability at Event Time
    # training
    Survival_Probs_event3_tr <- predict(Mod3, sel_clin_tr2,
                                        type ="prob_event")
    survivalProbs_p3_tr <- predict(Mod3, sel_clin_tr2,
                                   type = "prob_times")
    # test
    Survival_Probs_event3 <- predict(Mod3, sel_clin_te2, type = "prob_event")

    survivalProbs_p3 <- predict(Mod3, sel_clin_te2, type = "prob_times")

    ## create a data frame combining predicted mean, median, survival
    # probability and actual time and event
    surv_res3 <- cbind(meanSurv3, medianSurv3,
                       Survival_Probs_event3,
                       sel_clin_te2$OS_month, sel_clin_te2$OS)
    colnames(surv_res3) <- c("Mean", "Median_surv",
                             "Survival_Prob_event", "Actual_OS_time", "Event")

    # training data
    # extract prob times at diff times points
    survivalProbs_t_mat_3_tr <- as.matrix(survivalProbs_p3_tr)
    survivalProbs_t_mat1_3_tr <- survivalProbs_t_mat_3_tr[, -1]
    survivalProbs_t_mat1_t_3_tr <- t(survivalProbs_t_mat1_3_tr)
    survivalProbs_t_mat1_t2_3_tr <- survivalProbs_t_mat1_t_3_tr[, -1]

    # test
    # extract prob times at diff times points
    survivalProbs_t_mat_3 <- as.matrix(survivalProbs_p3)
    survivalProbs_t_mat1_3 <- survivalProbs_t_mat_3[, -1]
    survivalProbs_t_mat1_t_3 <- t(survivalProbs_t_mat1_3)
    survivalProbs_t_mat1_t2_3 <- survivalProbs_t_mat1_t_3[, -1]



    # Calcualte Evaluation parameters on training data
    # create survival object
    surv_obj_3_tr <- survival::Surv(sel_clin_tr2$OS_month, sel_clin_tr2$OS)

    # calculate IBS (Integrated Brier Score for test data
    IBS1_3_tr <- round(IBS(surv_obj_3_tr,
                           sp_matrix = survivalProbs_t_mat1_t2_3_tr,
                           survivalProbs_p3_tr$time[-1]), 3)

    # Concordance Index
    c_index_3_tr <- round(SurvMetrics::Cindex(surv_obj_3_tr,
                                              predicted = medianSurv3_tr), 2)

    # Combine evaluation parameters to get Matrix
    Error_mat_3_tr <- cbind(IBS1_3_tr, c_index_3_tr)
    # Calcualte Evaluation parameters on test data
    # create survival object
    surv_obj_3 <- survival::Surv(sel_clin_te2$OS_month,
                                 sel_clin_te2$OS)

    # calculate IBS (Integrated Brier Score for test data
    IBS1_3 <- round(IBS(surv_obj_3,
                        sp_matrix = survivalProbs_t_mat1_t2_3,
                        survivalProbs_p3$time[-1]), 3)

    # Concordance Index
    c_index_3 <- round(SurvMetrics::Cindex(surv_obj_3,
                                           predicted = medianSurv3), 2)
    # Combine evaluation parameters to get Matrix
    Error_mat_3_te <- cbind(IBS1_3, c_index_3)
    Error_mat_3 <- rbind(Error_mat_3_tr, Error_mat_3_te)
    #column names
    colnames(Error_mat_3) <- c("IBS_Score", "c_index")
    rownames(Error_mat_3) <- c("Training_set", "Test_set")
  } else if (Model_type == 4) { # model 4-  Univariate with Clin features
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
    # Prediction on test data
    # Survival curve of test data
    survCurves5 <- predict(Mod5, sel_clin_te2, type = "survivalcurve")
    # add column names
    colnames(survCurves5) <- c("time_point", rownames(sel_clin_te2))
    # make dataframe
    survCurves5_df <- as.data.frame(survCurves5)
    # Predicted Mean/Median Survival Time
    # Predicted Mean
    meanSurv5_tr <- predict(Mod5, sel_clin_tr2, type = "mean_time")
    meanSurv5 <- predict(Mod5, sel_clin_te2, type = "mean_time")
    # Median
    medianSurv5_tr <- predict(Mod5, sel_clin_tr2, type = "median_time")
    medianSurv5 <- predict(Mod5, sel_clin_te2, type = "median_time")
    # make data frames of predicted  mean and median survival time
    meanSurv5_d <- as.data.frame(meanSurv5)
    medianSurv5_d <- as.data.frame(medianSurv5)
    names_5 <- as.data.frame(rownames(sel_clin_te2))
    # combine predicted mean and median time
    mean_median_surv5_d <- cbind(names_5, meanSurv5_d, medianSurv5_d)
    # add row names and column names
    rownames(mean_median_surv5_d) <- c(rownames(sel_clin_te2))
    colnames(mean_median_surv5_d) <- c("IDs", "Mean", "Median")
    # Survival Probability at Event Time
    # training data set
    Survival_Probs_event5_tr <- predict(Mod5, sel_clin_tr2,
                                        type ="prob_event")
    survivalProbs_p5_tr <- predict(Mod5, sel_clin_tr2,
                                   type = "prob_times")
    # Test
    Survival_Probs_event5 <- predict(Mod5, sel_clin_te2,
                                     type = "prob_event")
    survivalProbs_p5 <- predict(Mod5, sel_clin_te2,
                                type = "prob_times")
    # create a data frame combining predicted mean, median, survival
    # probability and actual time and event
    surv_res5 <- cbind(meanSurv5, medianSurv5, Survival_Probs_event5,
                       sel_clin_te2$OS_month, sel_clin_te2$OS)
    # add column and row names
    colnames(surv_res5) <- c("Mean", "Median_surv",
                             "Survival_Prob_event",
                             "Actual_OS_time", "Event")
    rownames(surv_res5) <- rownames(sel_clin_te2)
    # extract prob times at diff times points on training data
    survivalProbs_t_mat_5_tr <- as.matrix(survivalProbs_p5_tr)
    survivalProbs_t_mat1_5_tr <- survivalProbs_t_mat_5_tr[, -1]
    survivalProbs_t_mat1_t_5_tr <- t(survivalProbs_t_mat1_5_tr)
    survivalProbs_t_mat1_t2_5_tr <- survivalProbs_t_mat1_t_5_tr[, -1]
    # extract prob times at diff times points on test data
    survivalProbs_t_mat_5 <- as.matrix(survivalProbs_p5)
    survivalProbs_t_mat1_5 <- survivalProbs_t_mat_5[, -1]
    survivalProbs_t_mat1_t_5 <- t(survivalProbs_t_mat1_5)
    survivalProbs_t_mat1_t2_5 <- survivalProbs_t_mat1_t_5[, -1]
    # Calcualte Evalulation parameters on Training data
    # create survival object
    surv_obj_5_tr <- survival::Surv(sel_clin_tr2$OS_month, sel_clin_tr2$OS)
    # calculate IBS (Integrated Brier Score for test data
    IBS1_5_tr <- round(IBS(surv_obj_5_tr,
                           sp_matrix = survivalProbs_t_mat1_t2_5_tr,
                           survivalProbs_p5_tr$time[-1]), 3)
    # Concordance Index
    c_index_5_tr <- round(SurvMetrics::Cindex(surv_obj_5_tr,
                                              predicted = medianSurv5_tr), 2)

    # Combine evaluation parameters to get Matrix
    Error_mat_5_tr <- cbind(IBS1_5_tr, c_index_5_tr)
    # Calcualte Evalulation parameters on test data
    # create survival object
    surv_obj_5 <- survival::Surv(sel_clin_te2$OS_month, sel_clin_te2$OS)
    # calculate IBS (Integrated Brier Score for test data
    IBS1_5 <- round(IBS(surv_obj_5,
                        sp_matrix = survivalProbs_t_mat1_t2_5,
                        survivalProbs_p5$time[-1]), 3)
    # Concordance Index
    c_index_5 <- round(SurvMetrics::Cindex(surv_obj_5,
                                           predicted = medianSurv5), 2)
    # Combine evaluation parameters to get Matrix
    Error_mat_5_te <- cbind(IBS1_5, c_index_5)
    Error_mat_5 <- rbind(Error_mat_5_tr, Error_mat_5_te)
    #add column and row names
    colnames(Error_mat_5) <- c("IBS_Score", "c_index")
    rownames(Error_mat_5) <- c("Training_set", "Test_set")
  }
  if (Model_type == 1) {
    survCurves_df <- survCurves1_df
    mean_median_surv_d <- mean_median_surv1_d
    Error_mat <- Error_mat_1
    surv_res <- surv_res1
  }
  if (Model_type == 2) {
    survCurves_df <- survCurves2_df
    mean_median_surv_d <- mean_median_surv2_d
    Error_mat <- Error_mat_2
    surv_res <- surv_res2
  }
  if (Model_type == 3) {
    survCurves_df <- survCurves3_df
    mean_median_surv_d <- mean_median_surv3_d
    Error_mat <- Error_mat_3
    surv_res <- surv_res3
  }
  if (Model_type == 4) {
    survCurves_df <- survCurves5_df
    mean_median_surv_d <- mean_median_surv5_d
    Error_mat <- Error_mat_5
    surv_res <- surv_res5
  }
  # Return a list containing data.
  return(list(
    survCurves_data = survCurves_df,
    mean_median_survival_time_data = mean_median_surv_d,
    survival_result_based_on_MTLR = surv_res,
    Error_mat_for_Model = Error_mat
  ))
}
