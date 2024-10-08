
#' This function selects significant features (p-value <0.05)s based on median
#' expression values of features using Univariate Survival analysis.
#' @param train_data :args1 - training data (Patients data with clinical and
#' gene expression, where samples are in rows and features/genes are in
#' columns)
#' @param test_data_data :args2 - training data (Patients data with clinical
#' and gene expression, where samples are in rows and features/genes are in
#' columns)
#' @param col_num :args3 - column number in data at where clinical info ends
#' @param surv_time :arg4 - name of column which contain survival time
#' (in days) information
#' @param surv_event :arg5 - name of column which contain survival
#' eventinformation
#' @param output_univariate_train :args6- name of output to store univariate
#' selected features for training data
#' @param output_univariate_test :args7- name of output to store univariate
#' selected features for test data
#' @import dplyr
#' @import survival
#' @import survminer
#' @import ggplot2
#' @examples
#' data(Train_Norm_data, package = "CPSM")
#' data(Test_Norm_data, package = "CPSM")
#' Univariate_sig_features_f(train_data=Train_Norm_data, test_data=
#' Test_Norm_data, col_num=21, surv_time="OS_month" , surv_event="OS")
#' Usage: Univariate_sig_features_f(train_data, test_data, col_num, surv_time,
#' surv_event)
#' @export



Univariate_sig_features_f <- function(train_data, test_data, col_num,surv_time,
                                      surv_event){
  # Check if any input variable is empty
  if (length(train_data) == 0 ||  length(test_data) == 0|| length(col_num) == 0
      ||length(surv_time) == 0 ||length(surv_event) == 0 ) {
    message("Error: Empty input variable detected.")
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

  #load data
  tr_data1 <- train_data
  te_data1 <- test_data

  # rename survival time and event column name
  colnames(tr_data1)[colnames(tr_data1) == surv_time] <- "OS_month"
  colnames(tr_data1)[colnames(tr_data1) == surv_event] <- "OS"

  # rename survival time and event column name
  colnames(te_data1)[colnames(te_data1) == surv_time] <- "OS_month"
  colnames(te_data1)[colnames(te_data1) == surv_event] <- "OS"

  ###################### Create Survival Object #####################
  surv_object <- Surv(time = tr_data1$OS_month, event = tr_data1$OS)
  # create a file to store results

  n <- col_num - 1
  ### Significant results
  #---------------------------------------------------------------------------
  # Initialize a list to store results
  results_list <- list()
  # Perform uni variate survival analysis for each feature based on median
  # expression
  for(i in seq(from = col_num, to = length(tr_data1), by = 1)) {
    # Create survival object
    surv_object <- Surv(time = tr_data1$OS_month, event = tr_data1$OS)

    # Survival analysis: fits cox ph model to find HR for median cut
    fit1 <- survfit(surv_object ~ (tr_data1[, i]) > (median(tr_data1[1, i])),
                    data = tr_data1)
    #summary(fit1)
    fit1.coxph <- coxph(surv_object ~ (tr_data1[, i]) >
                          (median(tr_data1[1, i])), data = tr_data1)
    # summary(fit1.coxph)
    first <- coef(summary(fit1.coxph))

    # Check whether the p-value is significant (< 0.05) or not
    if ((first[5] <= 0.05) && (!is.na(first[5])) && (!is.na(first[2]))) {

      # Store results in the list
      results_list[[length(results_list) + 1]] <- c(
        ID = colnames(tr_data1[i]),
        Beta = first[1],
        HR = first[2],
        `P-value` = first[5],
        GP1 = fit1$n[1],
        GP2 = fit1$n[2],
        `Hr-Inv-lst` = 1 / first[2],
        Concordance = fit1.coxph$concordance[6],
        Std_Error = fit1.coxph$concordance[7]
      )
    }
  }

  # Convert the list to a data frame for easier handling
  results_df <- do.call(rbind, results_list)

  # Set correct column names
  colnames(results_df) <- c("ID", "Beta", "HR", "P-value", "GP1", "GP2",
                            "Hr-Inv-lst", "Concordance", "Std_Error")
  selected_feature_names <- results_df[, 1]

  # Prepare training data with selected features
  sel_univ_train <- tr_data1[, colnames(tr_data1) %in% selected_feature_names,
                             drop = FALSE]
  # Convert to data frame if necessary
  sel_univ_train_1 <- as.data.frame(sel_univ_train)

  # Prepare test data with selected features
  sel_univ_test <- te_data1[, colnames(te_data1) %in% selected_feature_names,
                            drop = FALSE]
  # Convert to data frame if necessary
  sel_univ_test_1 <- as.data.frame(sel_univ_test)
  
  
  ### Significant results for clinical features
  #---------------------------------------------------------------------------
  # Initialize a list to store results
  results_list2 <- list()
  
  # Extract data for clinical features only
  tr_data_clin <-  tr_data1[1:n]
  te_data_clin <-  tr_data1[1:n]
  
  # Exclude OS and OS month columns by name
  tr_data2 <- tr_data_clin[, !colnames(tr_data_clin) %in% c("OS", "OS_month")]
  te_data2 <- te_data_clin[, !colnames(te_data_clin) %in% c("OS", "OS_month")]  
  
  # Perform uni variate survival analysis for each clinical feature 
  for(i in seq(from = 1, to = length(tr_data2), by = 1)) {    
    # Create survival object
    surv_object <- Surv(time = tr_data1$OS_month, event = tr_data1$OS)
    
    tryCatch({
      
      # Survival analysis: fits cox ph model to find HR for median cut
      fit2 <- survfit(surv_object ~ tr_data2[, i],data = tr_data2)
      #summary(fit1)
      fit2.coxph <- coxph(surv_object ~ tr_data2[, i], data = tr_data2)
      # summary(fit1.coxph)
      first2 <- coef(summary(fit2.coxph))
      
      # Check whether the p-value is significant (< 0.05) or not
      if ((first2[5] <= 0.05) && (!is.na(first2[5])) && (!is.na(first2[2]))) {
        
        # Store results in the list
        results_list2[[length(results_list2) + 1]] <- c(
          ID = colnames(tr_data2[i]),
          Beta = first2[1],
          HR = first2[2],
          `P-value` = first2[5],
          GP1 = fit2$n[1],
          GP2 = fit2$n[2],
          `Hr-Inv-lst` = 1 / first2[2],
          Concordance = fit2.coxph$concordance[6],
          Std_Error = fit2.coxph$concordance[7]
        )
      }
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n") })
  }
  
  # Convert the list to a data frame for easier handling
  results_df2 <- do.call(rbind, results_list2)
  
  # Set correct column names
  colnames(results_df2) <- c("ID", "Beta", "HR", "P-value", "GP1", "GP2",
                             "Hr-Inv-lst", "Concordance", "Std_Error")
  selected_feature_names2 <- results_df2[, 1]
  
  # Prepare training data with selected features
  sel_univ_train2 <- tr_data1[, colnames(tr_data1) %in% selected_feature_names2,
                              drop = FALSE]
  # Convert to data frame if necessary
  sel_univ_train_2 <- as.data.frame(sel_univ_train2)
  
  # Prepare test data with selected features
  sel_univ_test2 <- te_data2[, colnames(te_data2) %in% selected_feature_names2,
                             drop = FALSE]
  # Convert to data frame if necessary
  sel_univ_test_2 <- as.data.frame(sel_univ_test2)
  
  
  
  

  # Return a list containing data.
  return(list(Univariate_Survival_Significant_genes_List = results_df,
              Train_Uni_sig_data = sel_univ_train_1,
              Test_Uni_sig_data = sel_univ_test_1,
              Univariate_Survival_Significant_clin_List = results_df2,
              Train_Uni_sig_clin_data = sel_univ_train_2,
              Test_Uni_sig_clin_data = sel_univ_test_2 ))

}

