#' Nomogram Generation Function
#'
#' This function generates a Cox proportional hazards nomogram for survival
#' prediction based on user-defined features and survival information. It also
#' computes bias-corrected and original C-index values for model validation.
#'
#' @param data A data frame containing survival information and predictor
#' variables. Must include columns for survival time, event status, and
#' selected features.
#' @param Feature_List A data frame or list containing the selected features
#' (column names) to be included in the nomogram.
#' @param surv_time A character string specifying the column name for survival
#'  time in the input data.
#' @param surv_event A character string specifying the column name for survival
#' event status in the input data (1 for event, 0 for censored).
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{C_index_mat}: A matrix containing the bias-corrected and
#'   original C-index values.
#' }
#'
#' @details The function dynamically creates a Cox proportional hazards model
#' based on the selected features and survival data. It computes survival
#' probabilities for 1, 3, 5, and 10 years (if applicable) and plots the
#' nomogram. The function uses internal validation with bootstrapping
#' (B = 1000) to calculate the bias-corrected C-index.
#'
#' The following steps are performed:
#' \itemize{
#'   \item Validates input data and renames survival columns to standard names.
#'   \item Prepares the data for modeling using \code{rms::datadist}.
#'   \item Dynamically creates a Cox proportional hazards model formula.
#'   \item Fits the Cox model using \code{rms::cph}.
#'   \item Generates the nomogram with survival probabilities at specified
#'   time points.
#'   \item Validates the model and computes C-index values.
#' }
#'
#' @import reshape2
#' @import ggplot2
#' @import rms
#' @import survivalROC
#' @import Matrix
#' @import svglite
#'
#' @examples
#' data(Train_Data_Nomogram_input, package = "CPSM")
#' data(feature_list_for_Nomogram, package = "CPSM")
#' Nomogram_generate_f(
#'   data = Train_Data_Nomogram_input, Feature_List =
#'     feature_list_for_Nomogram, surv_time = "OS_month", surv_event = "OS"
#' )
#'
#' @export


Nomogram_generate_f <- function(data, Feature_List, surv_time, surv_event) {
  # Check if any input variable is empty
  if (length(data) == 0 || length(Feature_List) == 0 || length(surv_time)
  == 0 || length(surv_event) == 0) {
    message("Error: Empty input variable detected.")
  }

  # load data
  data <- data

  # rename survival time and event column name
  colnames(data)[colnames(data) == surv_time] <- "OS_month"
  colnames(data)[colnames(data) == surv_event] <- "OS"

  # Load user defined a list of features for which user want to develop
  # nomogram
  ftr_list <- Feature_List

  # create data frame with selected features (user provided list)
  sel_data <- as.data.frame(data[, colnames(data) %in% c(ftr_list$ID), ])

  # add survival information
  sel_data1 <- cbind(data["OS"], data["OS_month"], sel_data)

  # computes statistical summaries of predictors using datadist, it convert
  # data into format to develop model for nomogram.
  d <- sel_data1
  ddist <- rms::datadist(d)
  options(datadist = ddist)

  # Create the formula dynamically
  formula_str <- paste(
    "Surv(OS_month, OS == 1) ~",
    paste(colnames(d)[3:ncol(d)],
      collapse = " + "
    )
  )

  # Convert the character string to a formula
  formula_obj <- as.formula(formula_str)

  # Fit the cox nomogram model hazards model
  cox1 <- cph(formula_obj,
    x = TRUE, y = TRUE, data = d, surv = TRUE,
    time.inc = 3
  )

  ##### #################
  surv <- Survival(cox1)
  risk <- function(x) 1 / (1 + exp(-x))
  surv_1 <- function(x) surv(1 * 12, lp = x) # defined time.inc,1 year OS
  surv_2 <- function(x) surv(1 * 36, lp = x) # defined time.inc,3 year OS
  surv_3 <- function(x) surv(1 * 60, lp = x) # defined time.inc,5 year OS
  surv_4 <- function(x) surv(1 * 120, lp = x) # defined time.inc,10 year OS

  if (cox1$maxtime > 120) {
    nom_cox1 <- nomogram(cox1,
      fun = list(risk, surv_1, surv_2, surv_3, surv_4),
      lp = FALSE,
      funlabel = c(
        "Risk", "1-Year Survival Probability",
        "3-Year Survival Probability",
        "5-Year Survival Probability",
        "10-Year Survival Probability"
      ),
      maxscale = 100,
      fun.at = c(
        "1.0", "0.95", "0.90", "0.85", "0.80", "0.70",
        "0.6", "0.5", "0.4", "0.3", "0.2", "0.1"
      )
    )
  } else if (cox1$maxtime > 60 & cox1$maxtime < 120) {
    nom_cox1 <- nomogram(cox1,
      fun = list(risk, surv_1, surv_2, surv_3),
      lp = FALSE,
      funlabel = c(
        "Risk", "1-Year Survival Probability",
        "3-Year Survival Probability",
        "5-Year Survival Probability"
      ),
      maxscale = 100,
      fun.at = c(
        "1.0", "0.95", "0.90", "0.85", "0.80", "0.70",
        "0.6", "0.5", "0.4", "0.3", "0.2", "0.1"
      )
    )
  } else if (cox1$maxtime < 60) {
    nom_cox1 <- nomogram(cox1,
      fun = list(risk, surv_1, surv_2),
      lp = FALSE,
      funlabel = c(
        "Risk", "1-Year Survival Probability",
        "3-Year Survival Probability"
      ),
      maxscale = 100,
      fun.at = c(
        "1.0", "0.95", "0.90", "0.85", "0.80", "0.70",
        "0.6", "0.5", "0.4", "0.3", "0.2", "0.1"
      )
    )
  } else {
    nom_cox1 <- nomogram(cox1,
      fun = list(risk, surv_1, surv_2),
      lp = FALSE,
      funlabel = c(
        "Risk", "1-Year Survival Probability",
        "3-Year Survival Probability"
      ),
      maxscale = 100,
      fun.at = c(
        "1.0", "0.95", "0.90", "0.85", "0.80", "0.70",
        "0.6", "0.5", "0.4", "0.3", "0.2", "0.1"
      )
    )
  }

  plot(nom_cox1,
    xfrac = .2,
    font.size = 0.5,
    cex.axis = 0.5,
    force.label = TRUE,
    tcl = 0.5,
    lmgp = 0.5,
    vnames = "labels",
    col.grid = gray(c(0.85, 0.95))
  )

  # Redirect console output to a temporary file
  sink(tempfile())
  # Perform validation while suppressing output
  v <- suppressWarnings(suppressMessages(validate(cox1,
    dxy = TRUE, B = 1000
  )))
  # Restore console output
  sink()
  # Calculate C-index
  Dxy <- v[rownames(v) == "Dxy", colnames(v) == "index.corrected"]
  orig_Dxy <- v[rownames(v) == "Dxy", colnames(v) == "index.orig"]
  bias_corrected_c_index <- round(abs(Dxy) / 2 + 0.5, 2)
  orig_c_index <- round(abs(orig_Dxy) / 2 + 0.5, 2)

  C_index_mat <- cbind(bias_corrected_c_index, orig_c_index)
  colnames(C_index_mat) <- c("Bias-corrected C-index", "C-index")

  # Return a list containing data.
  return(list(C_index_mat = C_index_mat))
}
