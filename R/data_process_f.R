#' This function converting OS time (in days) into months and removing samples
#' where OS time information is missing. Note: In the Example data OS time
#' labelled is "OS.time" and OS event labelled as "OS"
#' @param data:args1 - data (Patients data with clinical and gene expression,
#' where samples are in rows and features/genes are in columns)
#' @param col_num:args2 - column number in data at where clinical info ends
#' @param surv_time:arg3 - name of column which contain survival time (in days)
#' information
#' @param output:arg4 - name of output file, in chich user want to store data
#' @import MASS
#' @import dplyr
#' @examples
#' data(Example_TCGA_LGG_FPKM_data, package = "CPSM")
#' data_process_f(
#'   data = assays(Example_TCGA_LGG_FPKM_data)$expression,
#'   col_num = 20, surv_time = "OS.time"
#' )
#' Usage:data_process_f(data, col_num, surv_time)
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
