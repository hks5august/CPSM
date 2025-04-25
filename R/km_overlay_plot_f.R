#' @title Overlay Kaplan-Meier Plot Function
#'
#' @description
#' This function generates an overlayed Kaplan-Meier (KM) survival plot to visualize the survival curves
#' of training data stratified by predicted risk groups (e.g., High Risk vs Low Risk), along with the
#' survival curve of a selected test sample. It helps to visually assess how the test sample's predicted
#' survival aligns with the survival profiles of the training groups.
#'
#' @name km_overlay_plot_f
#'
#' @param Train_results A data frame containing the training samples' survival prediction results (obtained from RF risk prediction model).
#'   The following columns are required:
#'   - `Sample_ID`: Unique identifier for the sample.
#'   - `Actual`: Actual TCGA ID (often the same as `Sample_ID`).
#'   - `Predicted_Risk_Group`: Predicted class label (e.g., `High_Risk` / `Low_Risk`).
#'   - `High_Risk_Prob`: Probability of being classified as High Risk.
#'   - `Low_Risk_Prob`: Probability of being classified as Low Risk.
#'   - `Prediction_Prob`: Maximum class probability (confidence in the prediction).
#'   - `OS_month`: Overall survival time in months.
#'   - `OS_event`: Event indicator (`1` = death, `0` = censored).
#'   The data frame should have no duplicate `Sample_ID`s, and the `OS_month` column should be numeric while `OS_event` should be an integer.
#'
#' @param Test_results A data frame containing the test samples' survival prediction results (obtained from RF risk prediction model).
#'   The structure is similar to `Train_results`, with the following columns:
#'   - `Sample_ID`: Unique identifier for the sample.
#'   - `Actual`: Actual TCGA ID (often the same as `Sample_ID`).
#'   - `Predicted_Risk_Group`: Predicted class label (e.g., `High_Risk` / `Low_Risk`).
#'   - `High_Risk_Prob`: Probability of being classified as High Risk.
#'   - `Low_Risk_Prob`: Probability of being classified as Low Risk.
#'   - `Prediction_Prob`: Maximum class probability (confidence in the prediction).
#'   - `OS_month`: Overall survival time in months.
#'   - `OS_event`: Event indicator (`1` = death, `0` = censored).
#'   The data frame should have no duplicate `Sample_ID`s, and the `OS_month` column should be numeric while `OS_event` should be an integer.
#'   The `Sample_ID` of the `selected_sample` must be present in this data frame.
#'
#' @param survcurve_te_data A data frame containing survival probability curves for test samples. It should include:
#'   - `time_point`: Time points for survival estimates.
#'   - `<Sample_ID columns>`: Survival probabilities for each sample at each time point.
#'
#' @param selected_sample The `Sample_ID` of the selected test sample for which the user wants to plot the overlayed Kaplan-Meier plot.
#'   The `Sample_ID` should be present in the `Test_results` data frame.
#'
#' @return A `ggplot` object (from `ggsurvplot`) showing:
#'   - Kaplan-Meier survival curves for the training dataset, stratified by predicted risk group.
#'   - An overlaid survival curve for the selected test sample (in a dashed green line).
#'   - Annotations displaying the test sample's predicted risk group and prediction probability.
#'
#' @importFrom survival Surv survfit
#' @importFrom survminer ggsurvplot
#' @importFrom ggplot2 geom_step annotate
#' @importFrom reshape2 melt
#'
#' @examples
#' # Example usage of the km_overlay_plot_f function
#'
#' # Example Train_results data frame
#' Train_results <- data.frame(
#'   Sample_ID = c("TCGA-TQ-A7RQ-01", "TCGA-TQ-A7RQ-02", "TCGA-TQ-A7RQ-03"),
#'   Actual = c("TCGA-TQ-A7RQ-01", "TCGA-TQ-A7RQ-02", "TCGA-TQ-A7RQ-03"),
#'   Predicted_Risk_Group = c("High_Risk", "Low_Risk", "High_Risk"),
#'   High_Risk_Prob = c(0.8, 0.3, 0.85),
#'   Low_Risk_Prob = c(0.2, 0.7, 0.15),
#'   Prediction_Prob = c(0.8, 0.7, 0.85),
#'   OS_month = c(12, 24, 6),
#'   OS_event = c(1, 0, 1)
#' )
#'
#' # Example Test_results data frame
#' Test_results <- data.frame(
#'   Sample_ID = c("TCGA-TQ-A7RQ-04", "TCGA-TQ-A7RQ-05"),
#'   Actual = c("TCGA-TQ-A7RQ-04", "TCGA-TQ-A7RQ-05"),
#'   Predicted_Risk_Group = c("Low_Risk", "High_Risk"),
#'   High_Risk_Prob = c(0.25, 0.9),
#'   Low_Risk_Prob = c(0.75, 0.1),
#'   Prediction_Prob = c(0.75, 0.9),
#'   OS_month = c(18, 30),
#'   OS_event = c(0, 1)
#' )
#'
#' # Example survival probability data frame for test samples
#' survCurves_data <- data.frame(
#'   time_point = c(0, 6, 12, 18, 24, 30),
#'   TCGA_TQ_A7RQ_04 = c(1, 0.85, 0.7, 0.5, 0.3, 0.1),
#'   TCGA_TQ_A7RQ_05 = c(1, 0.95, 0.85, 0.65, 0.45, 0.25)
#' )
#'
#' # Generate Kaplan-Meier plot with overlay
#' KM_plot_results <- km_overlay_plot_f(
#'   Train_results = Train_results,
#'   Test_results = Test_results,
#'   survcurve_te_data = survCurves_data,
#'   selected_sample = "TCGA-TQ-A7RQ-04"
#' )
#'
#'
#' @export


utils::globalVariables(c("time", "surv"))

km_overlay_plot_f <- function(Train_results, Test_results, survcurve_te_data, selected_sample) {

  # Check Train_results
  if (missing(Train_results) || is.null(Train_results) || !is.data.frame(Train_results)) {
    message("Error: 'Train_results' must be a non-null data frame.")
    return(NULL)
  }

  # Check Test_results
  if (missing(Test_results) || is.null(Test_results) || !is.data.frame(Test_results)) {
    message("Error: 'Test_results' must be a non-null data frame.")
    return(NULL)
  }

  # Check survcurve_te_data
  if (missing(survcurve_te_data) || is.null(survcurve_te_data) || !is.data.frame(survcurve_te_data)) {
    message("Error: 'survcurve_te_data' must be a non-null data frame.")
    return(NULL)
  }

  # Check selected_sample
  if (missing(selected_sample) || is.null(selected_sample) || selected_sample == "") {
    message("Error: 'selected_sample' is missing, NULL, or an empty string.")
    return(NULL)
  }

  # Check if selected_sample exists in Test_results
  if (!(selected_sample %in% rownames(Test_results))) {
    message("Error: selected_sample '", selected_sample, "' is NOT present in the Test_results data.")
    return(NULL)
  } else {
    message("selected_sample '", selected_sample, "' is present in the Test_results data.")
  }


  # Create survival object for Train data
  Train_results$surv_Tr <- survival::Surv(time = Train_results$OS_month, event = Train_results$OS_event)
  km_Train_results_fit <- survfit(surv_Tr ~ Predicted_Risk_Group, data = Train_results)

  # KM plot for Train data
  km_Train_results_plot <- ggsurvplot(
    km_Train_results_fit,
    data = Train_results,
    risk.table = TRUE,
    break.time.by = 6,
    pval = FALSE,
    censor = TRUE,
    palette = c("red", "blue"),
    ggtheme = theme_minimal(),
    surv.median.line = "hv",
    title = "Kaplan-Meier Plot: Prediction of Test Sample in Comparison to Risk Groups of Training Data"
  )

  # Reshape Test survival curve data
  Test_results_long <- reshape2::melt(
    survcurve_te_data,
    id.vars = "time_point",
    variable.name = "Patient",
    value.name = "Survival_Probability"
  )

  # Extract the selected test sample data
  test_sample <- Test_results_long[Test_results_long$Patient == selected_sample, ]

  # Get predicted risk group and probability
  selected_test_pred_risk <- Test_results$Predicted_Risk_Group[rownames(Test_results) == selected_sample]
  selected_test_pred_Prob <- Test_results$Prediction_Prob[rownames(Test_results) == selected_sample]

  # Create survival object for the test sample
  test_sample$event <- ifelse(test_sample$Survival_Probability < 1, 1, 0)
  test_sample_surv <- survival::Surv(time = test_sample$time_point, event = test_sample$event)
  km_test_fit <- survfit(test_sample_surv ~ 1)

  # Prepare survival curve data
  test_curve <- data.frame(
    time = km_test_fit$time,
    surv = km_test_fit$surv
  )

  # Calculate median survival time for test sample
  median_test_surv_time <- tryCatch({
    test_curve$time[which.min(abs(test_curve$surv - 0.5))]
  }, error = function(e) { NA })

  # Add test sample curve and median line to KM plot
  km_Train_results_plot2 <- km_Train_results_plot$plot +
    geom_step(
      data = test_curve,
      aes(x = time, y = surv),
      color = "darkgreen",
      size = 0.6,
      linetype = "dashed"
    ) +

    annotate(
      "text",
      x = max(test_curve$time) * 1.5,
      y = 0.8,
      label = paste("Sample:", selected_sample,
                    "\nPredicted Risk:", selected_test_pred_risk,
                    "\nPrediction Probability:", round(selected_test_pred_Prob, 3)),
      color = "darkgreen",
      fontface = "bold"
    )

  return(km_Train_results_plot2)
}

