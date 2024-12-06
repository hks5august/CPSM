test_that("multiplication works", {
  data(New_data, package = "CPSM")
  result <- tr_test_f(data = assays(New_data)$expression, fraction = 0.9)
  expect_s3_class(result$train_data, "data.frame")
})
