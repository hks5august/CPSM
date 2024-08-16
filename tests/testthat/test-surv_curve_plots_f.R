test_that("multiplication works", {
  data(survCurves_data, package = "SPM")
  p <- surv_curve_plots_f(Surv_curve_data = survCurves_data,
                     selected_sample = "TCGA-TQ-A8XE-01")
  expect_s3_class(p, "ggplot")
})
