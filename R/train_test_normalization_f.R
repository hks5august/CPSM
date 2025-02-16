#' Train-Test Data Normalization Function
#'
#' This function performs normalization on training and test datasets, which
#' include clinical data and gene expression values.
#' The function applies log transformation and quantile normalization to the
#' gene expression values.
#'
#' @param train_data A data frame containing the training dataset. The first
#' columns represent clinical data, and the rest represent gene expression
#' data.
#' @param test_data A data frame containing the test dataset. The structure
#' is the same as `train_data`.
#' @param col_num An integer specifying the starting column index of the gene
#' expression data in the datasets.
#'
#' @return A list containing four data frames:
#' \describe{
#'   \item{Train_Clin}{A data frame containing the clinical data from the
#'   training dataset.}
#'   \item{Test_Clin}{A data frame containing the clinical data from the
#'   test dataset.}
#'   \item{Train_Norm_data}{A data frame containing the combined clinical
#'   and normalized gene expression data for the training set.}
#'   \item{Test_Norm_data}{A data frame containing the combined clinical
#'   and normalized gene expression data for the test set.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Splits clinical and gene expression data based on `col_num`.
#'   \item Applies log2 transformation to the gene expression data.
#'   \item Performs quantile normalization on the log-transformed data
#'   for both training and test datasets.
#'   \item Combines the clinical data with the normalized gene expression
#'   data.
#' }
#'
#' @examples
#' # Example datasets
#' train_data <- data.frame(matrix(rnorm(100), nrow = 10))
#' test_data <- data.frame(matrix(rnorm(100), nrow = 10))
#' col_num <- 5
#' result <- train_test_normalization_f(train_data, test_data, col_num)
#' str(result$Train_Norm_data)
#' str(result$Test_Norm_data)
#'
#' @import preprocessCore
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
