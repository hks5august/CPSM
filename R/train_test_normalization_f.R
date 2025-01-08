#' This function first tranform FPKM data into log2-scale transformation,
#' followed by quantile normalization, where training data values used as
#' target matrix both for training and test data for quantile  normalization
#' @param train_data :args1 - training data (TSV) (Patients data with clinical
#' and gene expression, where samples are in rows and features/genes are in
#' columns)
#' @param test_data :args2 - test data (TSV) (Patients data with clinical and
#' gene expression, where samples are in rows and features/genes are in
#' columns)
#' @param col_num :args3 - column number in data at where clinical info ends
#' @param train_clin_data: args4- name of training data output which stores
#' only clinical information
#' @param test_clin_data: args5 - name of test data output which stores only
#' clinical information
#' @param train_Normalized_data_clin_data: args6 -  output file name - outfile
#' file containing training clinical data and normalized data  (which is
#' log-scaled followed by quantile normalized data).
#' @param test_Normalized_data_clin_data: args7- output filename - outfile file
#' containing test clinical data and normalized data  (which is log-scaled
#' followed by quantile normalized data).
#' @import MASS
#' @import dplyr
#' @import  preprocessCore
#' @examples
#' data(train_FPKM, package = "CPSM")
#' data(test_FPKM, package = "CPSM")
#' train_test_normalization_f(
#'   train_data = train_FPKM, test_data = test_FPKM,
#'   col_num = 21
#' )
#' Usage:train_test_normalization_f(train_data, test_data, col_num)

#' @export


train_test_normalization_f <- function(train_data, test_data, col_num) {
  # Check if any input variable is empty
  if (length(train_data) == 0 || length(test_data) == 0 ||
    length(col_num) == 0) {
    message("Error: Empty input variable detected.")
  }
  # Check if col_num is valid
  if (is.na(col_num) || col_num <= 0) {
    message("Error: col_num is missing or invalid.")
  }

  # training and test data contains clinical data in first 21 columns and rest
  # columns repression gene expression values
  training <- train_data
  testing <- test_data

  n <- col_num - 1
  # extract clinical data
  tr_clin <- training[seq_len(n)]
  te_clin <- testing[seq_len(n)]


  ## Extract Expression dat
  tr_exp <- training[col_num:ncol(training)]
  te_exp <- testing[col_num:ncol(testing)]

  ## log scale transformation
  tr_exp1 <- (tr_exp + 1)
  tr_log_mat <- round(log2(tr_exp1), 3)

  te_exp1 <- (te_exp + 1)
  te_log_mat <- round(log2(te_exp1), 3)

  # quantile normalization
  # samples in columns and genes in the rows

  # transpose train data
  ref <- as.data.frame(t(tr_log_mat))

  # transpose test data
  test <- as.data.frame(t(te_log_mat))

  # convert into matrix
  ref1 <- as.matrix(ref)
  test1 <- as.matrix(test)

  # normalize train data
  norm_ref <- round(normalize.quantiles(ref1), 3)

  # target
  target <- normalize.quantiles.determine.target(norm_ref)

  # qunatile normlize test data
  tt <- round(normalize.quantiles.use.target(test1, target), 3)

  # transpose quantile train data
  norm_ref_t <- as.data.frame(t(norm_ref))
  colnames(norm_ref_t) <- colnames(tr_exp)
  rownames(norm_ref_t) <- rownames(tr_exp)

  # transpose quantile test data
  test_t <- as.data.frame(t(tt))
  colnames(test_t) <- colnames(te_exp)
  rownames(test_t) <- rownames(te_exp)

  # combine clin and normalzed data
  Train_norm_data <- cbind(tr_clin, norm_ref_t)
  Test_norm_data <- cbind(te_clin, test_t)

  Train_Clin <- cbind(tr_clin)
  Test_Clin <- cbind(te_clin)
  Train_Norm_data <- cbind(Train_norm_data)
  Test_Norm_data <- cbind(Test_norm_data)

  # Return a list containing data.
  return(list(
    Train_Clin = Train_Clin, Test_Clin = Test_Clin,
    Train_Norm_data = Train_Norm_data, Test_Norm_data = Test_Norm_data
  ))
}
