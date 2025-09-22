
test_that("Assess mean_median_surv_barplot_f", {
  # Load the data from the CPSM package
  data(mean_median_survival_time_data, package = "CPSM")

  # Run the function with the selected sample
  p1 <- mean_median_surv_barplot_f(surv_mean_med_data = mean_median_survival_time_data,
                                   selected_sample = "TCGA-TQ-A8XE-01",
				   font_size = 12, font_color = "darkblue")

  # Check that the output is a list
  expect_type(p1, "list")

  # Ensure the list contains two ggplot objects
  expect_s3_class(p1$mean_med_all_pat, "ggplot")
  expect_s3_class(p1$highlighted_selected_pat, "ggplot")
})
