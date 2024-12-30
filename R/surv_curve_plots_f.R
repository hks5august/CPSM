#' This function generate individual survival curve for patients. Moreover,
#' user can also highlight one specific sample by providing sample IDs
#' @param Surv_curve_data :args1 - Survival probabilty data for survival curve
#' of patients
#' @param selected_sample:args2 - ID of a specific sample user wants to
#' highlights on the plot against the background of curves of other samples
#' @import MASS
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import svglite
#' @examples
#' data(survCurves_data, package = "CPSM")
#' surv_curve_plots_f(
#'   Surv_curve_data = survCurves_data, selected_sample =
#'     "TCGA-TQ-A8XE-01"
#' )
#' usgae:surv_curve_plots_f(Surv_curve_data, selected_sample)
#' @export

utils::globalVariables(c("Time", "Value", "Patient"))


surv_curve_plots_f <- function(Surv_curve_data, selected_sample) {
  # load data
  survCurves_data <- Surv_curve_data

  # reshape the data in the form of long matrix
  survCurves_m <- melt(survCurves_data, id = "time_point")

  # add column names
  colnames(survCurves_m) <- c("Time", "Patient", "Value")

  # create survival curve plot for all samples
  Surv_curv_plot_all_pat <- ggplot(survCurves_m, aes(
    x = Time, y = Value,
    group = Patient,
    color = Patient
  )) +
    geom_line() +
    labs(x = "Time in Months", y = "Survival Probability") +
    ggtitle("Survival Curves for Patients") +
    geom_hline(yintercept = 0.5, color = "black", linetype = "dashed") +
    # add dashed line corresonds to 0.5 probability
    theme(
      legend.position = "bottom", legend.box = "vertical",
      legend.title = element_text(), legend.key = element_blank(),
      legend.key.size = unit(3, "mm"), legend.text = element_text(size = 4)
    ) +
    guides(color = guide_legend(title = "Patients"))

  # Plot
  # print(Surv_curv_plot_all_pat)


  # selected patient ID
  Selected_patient <- selected_sample

  # Check if the sample is in the column names
  if (!(selected_sample %in% colnames(Surv_curve_data))) {
    message(
      "The sample", selected_sample,
      "is not present in the Test dataset.\n"
    )
  } else {
    message("The sample", selected_sample, "is present in the Test dataset.\n")
  }

  # create survival curve plot for selected patient
  Surv_curv_plot_all_pats_with_highlighting_one_pat <- ggplot(
    survCurves_m, aes(
      x = Time, y = Value, linetype = Patient,
      color = Selected_patient
    )
  ) +
    scale_y_continuous(breaks = seq(0, 1, 0.1)) +
    ggtitle("Survival Curves for Patients with Highlightited Patient") +
    geom_hline(yintercept = 0.5, color = "red", linetype = "dashed") +
    # Set linetype for all lines based on Patient
    geom_line(size = 0.5) +
    # Highlight one patient with a specific color
    scale_linetype_manual(values = rep(
      "solid",
      length(unique(survCurves_m$Patient))
    )) +
    scale_color_manual(values = c("black", "red")) +
    geom_line(aes(colour = "yellow"), size = 0.5, data = ~
      subset(survCurves_m, Patient == Selected_patient)) +
    labs(x = "Time in Months", y = "Survival Probability") +
    theme(
      legend.position = "bottom",
      legend.title = element_text(),
      legend.key = element_blank(),
      legend.key.size = unit(3, "mm"),
      legend.text = element_text(size = 4)
    ) +
    guides(
      linetype = guide_legend(title = "Patients"), color =
        guide_legend(title = "Selected Patient")
    )

  # Plot
  # print(Surv_curv_plot_all_pats_with_highlighting_one_pat )
  # Return the plots as a list
  return(list(
    all_patients_plot = Surv_curv_plot_all_pat,
    highlighted_patient_plot =
      Surv_curv_plot_all_pats_with_highlighting_one_pat
  ))
}
