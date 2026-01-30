#' Perform Univariate Survival Analysis
#'
#' @description
#' This function performs univariate survival analysis on the provided
#' datasets to identify significant features
#' (genes or clinical features) associated with survival outcomes. It applies
#' Cox proportional hazards models
#' to compute hazard ratios and p-values for median expression cut-offs or
#' clinical feature values.
#'
#' @param train_data A data frame containing the training dataset. It must
#' include survival time and event columns.
#' @param test_data A data frame containing the test dataset. It must include
#'  survival time and event columns.
#' @param col_num An integer specifying the column number in the dataset from
#' which feature analysis should begin.
#' @param surv_time A string specifying the name of the survival time column
#' in the dataset.
#' @param surv_event A string specifying the name of the survival event column
#' in the dataset.
#'
#' @return A list containing:
#' \item{Univariate_Survival_Significant_genes_List}{A data frame of
#' significant gene features with associated statistics.}
#' \item{Train_Uni_sig_data}{A data frame of selected significant gene
#' features for the training dataset.}
#' \item{Test_Uni_sig_data}{A data frame of selected significant gene
#' features for the test dataset.}
#' \item{Univariate_Survival_Significant_clin_List}{A data frame of
#' significant clinical features with associated statistics.}
#' \item{Train_Uni_sig_clin_data}{A data frame of selected significant
#' clinical features for the training dataset.}
#' \item{Test_Uni_sig_clin_data}{A data frame of selected significant
#'  clinical features for the test dataset.}
#' \item{ZPH_Genes}{A list of `cox.zph` objects for each gene feature,
#' containing Schoenfeld residuals used to assess proportional hazards assumption.}
#' \item{PH_Summary_Genes}{A data frame summarizing proportional hazards
#' test results for each gene, including p-values and indicators of
#' assumption violation.}
#'
#' @details
#' The function first checks the validity of the input variables, such as
#' whether the survival time and event columns
#' are present in the training and test datasets. It then renames these
#' columns to "OS_month" and "OS" internally
#' for consistency. For each feature starting from the specified column
#' number (`col_num`), univariate Cox proportional
#' hazards models are fitted to assess the association between the feature
#' and survival outcomes.
#'
#' For clinical features, a similar process is applied, but the analysis is
#' limited to the initial `col_num - 1` columns.
#'
#' Significant features are determined based on a p-value threshold of 0.05.
#' The function returns the significant
#' features and their corresponding hazard ratios, p-values, and other related
#' statistics for both training and test datasets.
#'
#' @examples
#' # Example data (requires real survival datasets for testing)
#' train_data <- data.frame(
#'   OS_month = c(10, 15, 20, 25),
#'   OS = c(1, 0, 1, 1),
#'   Gene1 = c(2.3, 3.4, 1.2, 2.5),
#'   Gene2 = c(0.5, 1.8, 0.9, 0.6),
#'   Clinical1 = c(1, 0, 1, 0)
#' )
#' test_data <- data.frame(
#'   OS_month = c(12, 18, 22, 28),
#'   OS = c(1, 1, 0, 1),
#'   Gene1 = c(2.0, 3.1, 1.0, 2.3),
#'   Gene2 = c(0.4, 1.5, 1.0, 0.7),
#'   Clinical1 = c(1, 1, 0, 0)
#' )
#' #result <- Univariate_sig_features_f(
#' #   train_data = train_data,
#' #   test_data = test_data,
#' #   col_num = 3,
#' #   surv_time = "OS_month",
#' #   surv_event = "OS"
#' #)
#'
#' @export



Univariate_sig_features_f <- function(train_data, test_data, col_num,
                                      surv_time, surv_event) {
  # Check if any input variable is empty
  if (length(train_data) == 0 || length(test_data) == 0 ||
    length(col_num) == 0 || length(surv_time) == 0 ||
    length(surv_event) == 0) {
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

  # load data
  tr_data1 <- train_data
  te_data1 <- test_data

  # rename survival time and event column name
  colnames(tr_data1)[colnames(tr_data1) == surv_time] <- "OS_month"
  colnames(tr_data1)[colnames(tr_data1) == surv_event] <- "OS"

  # rename survival time and event column name
  colnames(te_data1)[colnames(te_data1) == surv_time] <- "OS_month"
  colnames(te_data1)[colnames(te_data1) == surv_event] <- "OS"

  # Create Survival Object
  surv_object <- Surv(time = tr_data1$OS_month, event = tr_data1$OS)
  # create a file to store results

  n <- col_num - 1
  #--------- Significant results ------------#
  # Initialize a list to store results
  results_list <- list()
  # Initialize lists to store zph test results
  zph_results_genes <- list()
  # Perform uni variate survival analysis for each feature based on median
  for (i in seq(from = col_num, to = ncol(tr_data1), by = 1)) {
    
  # Skip features with no variation
  if (length(unique(tr_data1[, i][!is.na(tr_data1[, i])])) < 2) {
    zph_results_genes[[colnames(tr_data1)[i]]] <- NA
    next
  }

  # Create survival object
    surv_object <- Surv(time = tr_data1$OS_month, event = tr_data1$OS)

    # Survival analysis: fits cox ph model to find HR for median cut
    fit1 <- survfit(surv_object ~ (tr_data1[, i] > median(tr_data1[, i], na.rm = TRUE)),
      data = tr_data1
    )
    # fitcoxph model
    fit1.coxph <- coxph(surv_object ~ (tr_data1[, i] > median(tr_data1[, i], na.rm = TRUE)), data = tr_data1)
      
    first <- coef(summary(fit1.coxph))
      
     # Perform zph test and store results
   zph_results_genes[[colnames(tr_data1[i])]] <- tryCatch(
   suppressMessages(cox.zph(fit1.coxph)),
  error = function(e) NA
)
    # coeff

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
  

# zph_results_genes has been filled with cox.zph objects

# Extract PH test p-values safely
if(length(zph_results_genes) == 0) {
  zph_pvals_genes <- numeric(0)
} else {
  zph_pvals_genes <- sapply(zph_results_genes, function(x) {
if (inherits(x, "cox.zph") && "GLOBAL" %in% rownames(x$table))
  as.numeric(x$table["GLOBAL", "p"])
else NA_real_  
})
}

# Create PH summary data frame safely
if(length(zph_pvals_genes) > 0) {
  PH_Summary_Genes <- data.frame(
    Feature = names(zph_pvals_genes),
    PH_pvalue = zph_pvals_genes,
    PH_Violation = zph_pvals_genes < 0.05
  )
} else {
  PH_Summary_Genes <- data.frame(
    Feature = character(0),
    PH_pvalue = numeric(0),
    PH_Violation = logical(0)
  )
}

  # Convert the list to a data frame for easier handling
  results_df <- do.call(rbind, results_list)

  # Set correct column names
  colnames(results_df) <- c(
    "ID", "Beta", "HR", "P-value", "GP1", "GP2",
    "Hr-Inv-lst", "Concordance", "Std_Error"
  )

# Get feature names
selected_feature_names <- results_df[, 1]

# Exclude features violating PH assumption
if(nrow(PH_Summary_Genes) > 0) {
  selected_feature_names <- selected_feature_names[
    !PH_Summary_Genes$PH_Violation[
      match(selected_feature_names, PH_Summary_Genes$Feature)
    ]
  ]
}  
 
  # Prepare training data with selected features
  sel_univ_train <- tr_data1[, colnames(tr_data1) %in% selected_feature_names,
    drop = FALSE
  ]
  # Convert to data frame if necessary
  sel_univ_train_1 <- as.data.frame(sel_univ_train)

  # Prepare test data with selected features
  sel_univ_test <- te_data1[, colnames(te_data1) %in% selected_feature_names,
    drop = FALSE
  ]
  # Convert to data frame if necessary
  sel_univ_test_1 <- as.data.frame(sel_univ_test)


  #---- Significant results for clinical features -------#
  #---------------------------------------------------------------------------
  # Initialize a list to store results
  results_list2 <- list()

  # Extract data for clinical features only
  tr_data_clin <- tr_data1[seq_len(n)]
  te_data_clin <- te_data1[seq_len(n)]


  # Exclude OS and OS month columns by name
  tr_data2 <- tr_data_clin[, !colnames(tr_data_clin) %in% c("OS", "OS_month")]
  te_data2 <- te_data_clin[, !colnames(te_data_clin) %in% c("OS", "OS_month")]

  # Perform uni variate survival analysis for each clinical feature
  for (i in seq(from = 1, to = ncol(tr_data2), by = 1)) {
    
   # Skip features with no variation
  if (length(unique(tr_data2[, i][!is.na(tr_data2[, i])])) < 2) {
    zph_results_genes[[colnames(tr_data2)[i]]] <- NA
    next
  }

    # Create survival object
    surv_object2 <- Surv(time = tr_data2$OS_month, event = tr_data2$OS)

    tryCatch(
      {
        # Survival analysis: fits cox ph model to find HR for median cut
        fit2 <- survfit(surv_object2 ~ tr_data2[, i], data = tr_data2)
        # COXPH model
        fit2.coxph <- coxph(surv_object2 ~ tr_data2[, i], data = tr_data2)
        # Coeff
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
      },
      error = function(e) {
        message("Massage :", conditionMessage(e), "\n")
      }
    )
  }


  # Convert the list to a data frame for easier handling
  results_df2 <- do.call(rbind, results_list2)

  # Set correct column names
  colnames(results_df2) <- c(
    "ID", "Beta", "HR", "P-value", "GP1", "GP2",
    "Hr-Inv-lst", "Concordance", "Std_Error"
  )
  selected_feature_names2 <- results_df2[, 1]

  # Prepare training data with selected features
  sel_univ_train2 <- tr_data1[, colnames(tr_data1) %in%
    selected_feature_names2,
  drop = FALSE
  ]
  # Convert to data frame if necessary
  sel_univ_train_2 <- as.data.frame(sel_univ_train2)

  # Prepare test data with selected features
  sel_univ_test2 <- te_data2[, colnames(te_data2) %in% selected_feature_names2, drop = FALSE]
  # Convert to data frame if necessary
  sel_univ_test_2 <- as.data.frame(sel_univ_test2)

  # Return a list containing data.
  return(list(
    Univariate_Survival_Significant_genes_List = results_df,
    Train_Uni_sig_data = sel_univ_train_1,
    Test_Uni_sig_data = sel_univ_test_1,
    Univariate_Survival_Significant_clin_List = results_df2,
    Train_Uni_sig_clin_data = sel_univ_train_2,
    Test_Uni_sig_clin_data = sel_univ_test_2 ,
    ZPH_Genes = zph_results_genes,
    PH_Summary_Genes = PH_Summary_Genes
  ))
}
