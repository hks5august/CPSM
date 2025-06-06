#' @title Mean and Median Survival Bar Plot
#'
#' @description This function generate barplots for mean and median survival
#' of patients. Besides, user can also highlight one specific sample by
#' providing sample IDs
#'
#' @name mean_median_surv_barplot_f
#'
#' @param surv_mean_med_data :args1 - Predicted Mean median survival time data
#' for patients
#' @param selected_sample :args2 - ID of the sample for which user want to
#' highlight in the plot.
#'
#' @return A list containing two `ggplot` objects:
#' \itemize{
#'   \item \code{mean_med_all_pat}: Bar plot of mean and median survival
#'   times for
#'   all patients.
#'   \item \code{highlighted_selected_pat}: Bar plot highlighting the selected
#'   patient's
#'   mean and median survival times with distinct colors.
#' }
#'
#' @import MASS
#' @import reshape2
#' @import ggplot2
#' @examples
#' data(mean_median_survival_time_data, package = "CPSM")
#' mean_median_surv_barplot_f(
#'   surv_mean_med_data =
#'     mean_median_survival_time_data,
#'   selected_sample = "TCGA-TQ-A8XE-01"
#' )
#'
#' @export
utils::globalVariables(c("IDs", "value", "variable"))



# Mean Survival and Median Survival time Barplots
mean_median_surv_barplot_f <- function(surv_mean_med_data, selected_sample) {
  mean_median_surv_d <- surv_mean_med_data
  # reshape data
  mean_median_surv_d_m <- melt(mean_median_surv_d)
  Barplot_mean_med_all_pat_surv <- ggplot(
    mean_median_surv_d_m,
    aes(
      x = .data$IDs, y = .data$value,
      fill = .data$variable,
      colour = .data$variable
    )
  ) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_y_continuous(breaks = seq(0, 150, 12)) +
    ggtitle("Predicted Mean/Median Survival Time of Patients") +
    labs(x = "Patients", y = "Predicted Survival time in Months") +
    theme(axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1,
      size = 6
    ))

  # Highlight Selected patient
  Selected_patient <- selected_sample
  # Create the barplot for selected patient and other patients with diff colors
  Barplot_with_highlighted_selected_pat <- ggplot(
    mean_median_surv_d_m,
    aes(
      x = IDs, y = value,
      fill = interaction(IDs == Selected_patient, variable)
    )
  ) +
    geom_bar(stat = "identity", position = "dodge") +
    # Add text labels on top of bars for the selected patient
    geom_text(
      aes(label = ifelse(IDs == Selected_patient, round(value, 2),
        ""
      )),
      vjust = -0.1, hjust = -0.1, color = "black"
    ) +

    # Customize the fill colors with a legend
    scale_fill_manual(
      values = c(
        "TRUE.Mean" = "red", "TRUE.Median" = "blue",
        "FALSE.Mean" = "lightgray", "FALSE.Median" = "darkgray"
      ),
      labels = c(
        "TRUE.Mean" = "Selected Patient Mean",
        "TRUE.Median" = "Selected Patient Median",
        "FALSE.Mean" = "Other Patients Mean",
        "FALSE.Median" = "Other Patients Median"
      ),
      name = "Legend"
    ) +

    # Add labels and other formatting as needed
    labs(
      title = "Predicted Mean/Median Survival Time of Patients",
      x = "Patients",
      y = "Survival Time (Months)"
    ) +
    theme_minimal() +
    scale_y_continuous(breaks = seq(0, 160, 12)) +
    theme(axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1,
      size = 8
    ))

  # Display the plot
  # Return the plots as a list
  return(list(
    mean_med_all_pat = Barplot_mean_med_all_pat_surv,
    highlighted_selected_pat = Barplot_with_highlighted_selected_pat
  ))
}
