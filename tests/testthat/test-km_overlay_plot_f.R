test_that("multiplication works", {
data(Train_results, package = "CPSM")
data(Test_results, package = "CPSM")
data(survCurves_data, package = "CPSM")
p <-  km_overlay_plot_f(Train_results = Train_results, Test_results = Test_results, 
survcurve_te_data = survCurves_data,  selected_sample = "TCGA-TQ-A7RQ-01") 

# Check that the output is a list
  expect_type(p, "list")

  # Ensure the list contains a ggplot object
  expect_s3_class(p, "ggplot")
})



