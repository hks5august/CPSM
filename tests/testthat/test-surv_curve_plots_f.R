test_that("multiplication works", {
  data(survCurves_data, package = "CPSM")
  p <- surv_curve_plots_f(Surv_curve_data = survCurves_data,
                          selected_sample = "TCGA-TQ-A8XE-01")
  # Check that the output is a list
  expect_type(p, "list")

  # Ensure the list contains two ggplot objects
  expect_s3_class(p$all_patients_plot, "ggplot")
  expect_s3_class(p$highlighted_patient_plot, "ggplot")
})



