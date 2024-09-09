test_that("Assess Univariate_sig_features_f", {
  data(Train_Norm_data, package = "CPSM")
  data(Test_Norm_data, package = "CPSM")
  Result_Uni <- Univariate_sig_features_f(train_data = Train_Norm_data,
                                          test_data = Test_Norm_data,
                                          col_num=21,
                                          surv_time = "OS_month" ,
                                          surv_event = "OS")
  Univariate_Su <- Result_Uni$Train_Uni_sig_data
  expect_s3_class(Univariate_Su, "data.frame")
})
