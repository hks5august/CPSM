
#' This function split data into training and test data as per user defined
#' ratio
#' @param data:args1 a data frame where features in the columns and samples
#'  must be in rows
#' @param fraction:args2 , by which user want to split data for training data;
#' e.g. 90% training, so fraction would be 0.9
#' @param train_data:args3, name of training set that user can provide
#' @param test_data:args4, name of test set that user can provide
#' @return  training and test data
#' @import MASS
#' @import dplyr
#' @examples
#' data(New_data, package = "CPSM")
#' tr_test_f(data=New_data, fraction=0.9)
#' Usage: tr_test_f(data, fraction)
#' @export

tr_test_f <- function(data, fraction)
{
  # Check if any input variable is empty
  if (is.null(data) || nrow(data) == 0) {
    message("Error: Data is NULL or empty.")
  }
  # Check if fraction or filenames are missing
  if (is.na(fraction) || fraction <= 0 || fraction >= 1) {
    message("Error: Fraction must be a number between 0 and 1.")
  }
  #define train-test split fraction
  numberTrain <- floor(nrow(data)* fraction)
  #extract training index
  #trInd <- sample(1:nrow(data), numberTrain)
  trInd <- sample(seq_len(nrow(data)), numberTrain)
  #create training data
  training <- data[trInd,]
  #create test data
  testing <- data[-trInd,]
  ######################### Write into files  #################################
  train_data <- data.frame('ID'=rownames(training), training)
  test_data <- data.frame('ID'=rownames(testing), testing)
  # Return a list containing train_data and test_data
  return(list(train_data = train_data, test_data = test_data))
}






