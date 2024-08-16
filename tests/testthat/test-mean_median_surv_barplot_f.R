test_that("Assess mean_median_surv_barplot_f", {
  data(mean_median_survival_time_data, package = "SPM")
  p1 <- mean_median_surv_barplot_f(surv_mean_med_data = mean_median_survival_time_data,
                             selected_sample = "TCGA-TQ-A8XE-01")
  expect_s3_class(p1, "ggplot")
})
