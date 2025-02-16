#' Train-Test Split Function
#'
#' This function splits a given dataset into training and testing sets based
#' on a specified fraction.
#'
#' @param data A data frame containing the dataset to be split.
#' @param fraction A numeric value between 0 and 1 indicating the proportion
#' of data to include in the training set.
#'
#' @return A list containing two data frames:
#' \describe{
#'   \item{train_data}{A data frame containing the training set.}
#'   \item{test_data}{A data frame containing the testing set.}
#' }
#'
#' @details
#' The function checks whether the input `data` is valid and whether the
#' `fraction` is a number between 0 and 1.
#' It then performs a random split of the dataset into training and testing
#' sets.
#'
#' @examples
#' # Example dataset
#' data <- data.frame(x = rnorm(100), y = rnorm(100))
#' fraction <- 0.7
#' result <- tr_test_f(data, fraction)
#' head(result$train_data)
#' head(result$test_data)
#'
#' @export

tr_test_f <- function(data, fraction) {
  # Check if any input variable is empty
  if (is.null(data) || nrow(data) == 0) {
    message("Error: Data is NULL or empty.")
  }
  # Check if fraction or filenames are missing
  if (is.na(fraction) || fraction <= 0 || fraction >= 1) {
    message("Error: Fraction must be a number between 0 and 1.")
  }
  # define train-test split fraction
  numberTrain <- floor(nrow(data) * fraction)
  # extract training index
  trInd <- sample(seq_len(nrow(data)), numberTrain)

  # create training data
  training <- data[trInd, ]
  # create test data
  testing <- data[-trInd, ]
  # Write into files
  train_data <- data.frame(training)
  test_data <- data.frame(testing)
  # Return a list containing train_data and test_data
  return(list(train_data = train_data, test_data = test_data))
}
