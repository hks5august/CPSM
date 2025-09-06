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

test_that("km_overlay_plot_f works with custom aesthetics", {
  data(Train_results, package = "CPSM")
  data(Test_results, package = "CPSM")
  data(survCurves_data, package = "CPSM")

  p <- km_overlay_plot_f(
    Train_results = Train_results,
    Test_results = Test_results,
    survcurve_te_data = survCurves_data,
    selected_sample = "TCGA-TQ-A7RQ-01",
    font_size = 16,
    line_size = 1.2,
    train_palette = c("purple", "orange"),
    test_line_col = "black",
    test_line_type = "solid"
  )

  # Still returns a valid ggplot
  expect_s3_class(p, "ggplot")
})

