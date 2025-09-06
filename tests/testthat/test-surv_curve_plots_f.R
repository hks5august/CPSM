test_that("surv_curve_plots_f works with defaults", {
  data(survCurves_data, package = "CPSM")
  p <- surv_curve_plots_f(Surv_curve_data = survCurves_data,
                          selected_sample = "TCGA-TQ-A8XE-01")

  # Check that the output is a list
  expect_type(p, "list")

  # Ensure the list contains two ggplot objects
  expect_s3_class(p$all_patients_plot, "ggplot")
  expect_s3_class(p$highlighted_patient_plot, "ggplot")
})

test_that("surv_curve_plots_f applies custom styling arguments", {
  data(survCurves_data, package = "CPSM")
  p <- surv_curve_plots_f(
    Surv_curve_data = survCurves_data,
    selected_sample = "TCGA-TQ-A8XE-01",
    font_size = 16,
    line_size = 2,
    all_line_col = "blue",
    highlight_col = "orange"
  )

  # Still should return list of ggplot objects
  expect_type(p, "list")
  expect_s3_class(p$all_patients_plot, "ggplot")
  expect_s3_class(p$highlighted_patient_plot, "ggplot")

  # Check that custom color is applied to highlight
  # (pull data and check mapping)
  highlighted_layer <- p$highlighted_patient_plot$layers[[2]]
  expect_true("orange" %in% unlist(highlighted_layer$aes_params))
})


