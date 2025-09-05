test_that("multiplication works", {
  data(survCurves_data, package = "CPSM")

# Call the function
  p <- surv_curve_plots_f(
    Surv_curve_data = Surv_curve_data,
    selected_sample = "TCGA_TQ_A8XE_01",
    font_size = 12,
    line_size = 0.5,
    all_line_col = "black",
    highlight_col = "red"
  )

  # Check that the output is a list
  expect_type(p, "list")

  # Ensure the list contains two ggplot objects
  expect_s3_class(p$all_patients_plot, "ggplot")
  expect_s3_class(p$highlighted_patient_plot, "ggplot")
})



