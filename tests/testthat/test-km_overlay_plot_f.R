test_that("km_overlay_plot_f works with default parameters", {
  data(Train_results, package = "CPSM")
  data(Test_results, package = "CPSM")
  data(survCurves_data, package = "CPSM")

  p <- km_overlay_plot_f(
    Train_results = Train_results,
    Test_results = Test_results,
    survcurve_te_data = survCurves_data,
    selected_sample = "TCGA-TQ-A7RQ-01"
  )

  # Ensure the output is a ggplot object
  expect_s3_class(p, "ggplot")
})

test_that("km_overlay_plot_f works with custom styling parameters", {
  data(Train_results, package = "CPSM")
  data(Test_results, package = "CPSM")
  data(survCurves_data, package = "CPSM")

  p <- km_overlay_plot_f(
    Train_results     = Train_results,
    Test_results      = Test_results,
    survcurve_te_data = survCurves_data,
    selected_sample   = "TCGA-TQ-A7RQ-01",
    font_size         = 16,
    train_palette     = c("red", "blue"),
    test_curve_col    = "darkgreen",
    test_curve_size   = 2,
    test_curve_lty    = "dotted",
    annotation_col    = "green"
  )

  # Ensure output is still a ggplot object
  expect_s3_class(p, "ggplot")
})

