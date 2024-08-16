#' This function generate barplots for mean and median survival of patients.
#' Besides, user can also highlight one specific sample by providing sample IDs
#' @param surv_mean_med_data :args1 - Predicted Mean median survival time data
#' for patients
#' @param selected_sample :args2 - ID of the sample for which user want to
#' highlight in the plot
#' @import MASS
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import svglite
#' @examples
#' data(mean_median_survival_time_data, package = "CPSM")
#' mean_median_surv_barplot_f(surv_mean_med_data=
#' mean_median_survival_time_data,
#' selected_sample="TCGA-TQ-A8XE-01")
#' usgae: mean_median_surv_barplot_f(surv_mean_med_data, selected_sample)
#' @export
utils::globalVariables(c("IDs", "value", "variable"))



############ Mean Survival and Median Survival time Barplots ############
mean_median_surv_barplot_f <- function(surv_mean_med_data, selected_sample)
{
  mean_median_surv_d <- surv_mean_med_data
  #reshape data
  mean_median_surv_d_m <- melt(mean_median_surv_d )
  Barplot_mean_med_all_pat_surv <- ggplot(mean_median_surv_d_m,
                                          aes(x = .data$IDs, y = .data$value,
                                              fill = .data$variable,
                                              colour = .data$variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_y_continuous(breaks=seq(0,150,12)) +
    ggtitle("Predicted Mean/Median Survival Time of Patients") +
    labs(x = "Patients", y = "Predicted Survival time in Months") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                     size = 6))

  # print(Barplot_mean_med_all_pat_surv)
  ### Highlight Selected patient
  Selected_patient <- selected_sample

  # Create a bar plot with ggplot2
  Barplot_with_highlighted_selected_pat<- ggplot(mean_median_surv_d_m,
                                                 aes(x = IDs, y = value,
                                                     fill = variable,
                                                     colour="gray")) +
    geom_bar(aes(
      color = ifelse(IDs == Selected_patient, "Selected Sample",
                     "Other Samples")),
      linetype = ifelse(mean_median_surv_d_m$IDs == Selected_patient,
                        "dashed", "solid"),
      stat = "identity", position = "dodge") +

    # Add text labels on top of bars with conditional formatting
    geom_text(aes(label=ifelse(IDs == Selected_patient, round(value, 2), "")),
              vjust =  -0.1,  hjust = -0.1, color = "black") +
    # Manually specify outline color for the highlighted sample
    scale_color_manual(values = c( "white", "black")) +

    # Add labels or other formatting as needed
    labs(title = "Bar Plot with Highlighted Sample") +
    theme_minimal() +labs(x="Patients",y="Predicted Survival time in Months")+
    scale_y_continuous(breaks=seq(0,160,12)) +
    ggtitle("Predicted Mean/Median Survival Time of Patients")  +

    # Manually specify linetype for the highlighted sample
    scale_linetype_manual(name ="Outline Type", values= c("dashed", "solid"))+
    # Hide legend for linetype based on the highlighted sample
    theme(legend.key = element_blank(), legend.title = element_blank(),
          legend.position = "top", legend.box = "horizontal") +
    theme(axis.text.x = element_text(angle = 90, vjust= 0.5, hjust=1, size=8))

  # Plot
  # print(Barplot_with_highlighted_selected_pat)
  # Return the plots as a list
  return(list(
    mean_med_all_pat = Barplot_mean_med_all_pat_surv,
    highlighted_selected_pat = Barplot_with_highlighted_selected_pat
  ))

}


