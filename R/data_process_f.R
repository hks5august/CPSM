#' @title Data Processing f
#'
#' @description This function converting OS time (in days) into months and
#' removing samples where OS time information is missing. Note: In the Example
#' data OS time labelled is "OS.time" and OS event labelled as "OS"
#'
#' @details This function converting OS time (in days) into months and
#' removing samples where OS time information is missing. Note: In the Example
#' data OS time labelled is "OS.time" and OS event labelled as "OS"
#'
#' @name data_process_f
#'
#' @param data A data frame containing clinical and expression data.
#' @param col_num An integer indicating the column number where expression data
#' starts.
#' @param surv_time A character string specifying the name of the survival time
#' column in the data.
#' @return A data frame combining processed clinical and expression data.
#' Includes a new column `OS_month` for survival time in months, and retains
#' only rows with positive `OS_month`.
#' @import MASS
#' @import SummarizedExperiment
#' @examples
#' library(SummarizedExperiment)
#' data(Example_TCGA_LGG_FPKM_data, package = "CPSM")
#' data_process_f(
#'   data = assays(Example_TCGA_LGG_FPKM_data)$expression,
#'   col_num = 20, surv_time = "OS.time"
#' )
#'
#' @export
utils::globalVariables(c("OS_month"))


data_process_f <- function(data, col_num, surv_time) {
  # Check if any input variable is empty or missing
  if (is.null(data) || nrow(data) == 0 || is.null(col_num) ||
    is.null(surv_time)) {
    message("Error: Input variable is null or data frame is empty.")
  }
  if (any(is.na(col_num)) || any(is.na(surv_time))) {
    message("Error: Missing values in col_num, surv_time.")
  }
  if (!surv_time %in% colnames(data)) {
    message("Error: surv_time column not found in the data.")
  }
  n <- col_num - 1
  # Extract Clinical data
  # data_clin <- data[1:n]
  data_clin <- data[seq_len(n)]

  OS.time1 <- surv_time
  # Access column using column name
  column_data <- data_clin[[OS.time1]]
  # Convert days into months (OS time)
  data_clin$OS_month <- round(column_data / 30.417, 0)
  ## Extract Expression data
  data_exp <- data[col_num:ncol(data)]
  # Combine clinical and Expression data
  data1 <- cbind(data_clin, data_exp)
  data2 <- subset(data1, OS_month > 0)
  New_data <- data.frame(data2)
  return(New_data)
}
