test_that("multiplication works", {
  data(New_data, package = "CPSM")
  result <- tr_test_f(data = New_data, fraction = 0.9)
  expect_s3_class(result$train_data, "data.frame")
})
