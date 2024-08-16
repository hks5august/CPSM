test_that("Assess train_test_normalization_f", {
  data(train_FPKM, package = "SPM")
  data(test_FPKM, package = "SPM")
  Result_N_data <- train_test_normalization_f(train_data = train_FPKM,
                                              test_data = test_FPKM,
                                              col_num = 21)
   Train_Clin <- Result_N_data$Train_Clin
   expect_s3_class(Train_Clin, "data.frame")
})
