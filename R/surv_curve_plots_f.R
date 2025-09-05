#' Generate Survival Curve Plots
#'
#' This function creates survival curve plots for all patients and highlights a
#' specific patient based on the selected sample.
#'
#' @param Surv_curve_data A data frame containing survival curve data. The
#' data frame must include a column named `time_point` representing time and
#' other columns corresponding to patient survival probabilities.
#' @param selected_sample A character string specifying the patient ID to be
#' highlighted in the plot.
#'
#' @return A list containing two ggplot2 objects:
#'   \item{all_patients_plot}{A ggplot object displaying survival curves for
#'   all patients.}
#'   \item{highlighted_patient_plot}{A ggplot object displaying survival
#'   curves for all patients with the selected patient highlighted.}
#'
#' @details
#' - The input data (`Surv_curve_data`) must be structured such that the
#' first column (`time_point`) represents time points, and each subsequent
#' column represents the survival probabilities for a specific patient.
#' - If the `selected_sample` is not present in the column names
#' of `Surv_curve_data`, a message will indicate that the sample is missing.
#' - The function uses the `ggplot2` package for generating survival
#' curve plots.
#'
#' @import Hmisc
#' @import ggfortify
#'
#' @name surv_curve_plots_f
#'
#' @examples
#' # Example survival curve data
#' Surv_curve_data <- data.frame(
#'   time_point = c(0, 1, 2, 3, 4),
#'   Patient1 = c(1, 0.9, 0.8, 0.7, 0.6),
#'   Patient2 = c(1, 0.85, 0.75, 0.65, 0.55)
#' )
#'
#' # Generate plots with Patient1 highlighted
#' plots <- surv_curve_plots_f(Surv_curve_data, selected_sample = "Patient1", font_size = 12, line_size = 0.5, all_line_col = "black", highlight_col = "red")
#'
#' # View the plots
#' print(plots$all_patients_plot)
#' print(plots$highlighted_patient_plot)
#'
#' @import ggplot2
#' @import reshape2
#'
#' @export


utils::globalVariables(c("Time", "Value", "Patient"))


surv_curve_plots_f <- function(Surv_curve_data, selected_sample, font_size = 12, line_size = 0.5,
                               all_line_col = "black", highlight_col = "red") {
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
    geom_line(size = line_size) +
    labs(x = "Time in Months", y = "Survival Probability") +
    ggtitle("Survival Curves for Patients") +
    geom_hline(yintercept = 0.5, color = "black", linetype = "dashed") +
    # add dashed line corresonds to 0.5 probability
    theme(
      legend.position = "bottom", legend.box = "vertical",
      legend.title = element_text(size = font_size), legend.key = element_blank(),
      legend.key.size = unit(3, "mm"), legend.text = element_text(size = font_size)
    ) +
    guides(color = guide_legend(title = "Patients"))

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
    geom_line(size = line_size) +
    # Highlight one patient with a specific color
    scale_linetype_manual(values = rep(
      "solid",
      length(unique(survCurves_m$Patient))
    )) +
    scale_color_manual(values = c("black", "red")) +
    geom_line(aes(colour = highlight_col),
      size = line_size,
      data = ~ subset(
        survCurves_m,
        Patient == Selected_patient
      )
    ) +
    labs(x = "Time in Months", y = "Survival Probability") +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = font_size),
      legend.key = element_blank(),
      legend.key.size = unit(3, "mm"),
      legend.text = element_text(size = font_size)
    ) +
    guides(
      linetype = guide_legend(title = "Patients"),
      color = guide_legend(title = "Selected Patient")
    )

  # Return the plots as a list
  return(list(
    all_patients_plot = Surv_curv_plot_all_pat,
    highlighted_patient_plot =
      Surv_curv_plot_all_pats_with_highlighting_one_pat
  ))
}
