test_that("multiplication works", {
  data(Train_Norm_data, package = "CPSM")
  data(Test_Norm_data, package = "CPSM")
  Result_PI <- Lasso_PI_scores_f(train_data = Train_Norm_data,
                                 test_data = Test_Norm_data,
                                 nfolds=5,
                                 col_num=21,
                                 surv_time = "OS_month",
                                 surv_event = "OS")
  expect_s3_class(Result_PI$Train_Lasso_key_variables, "data.frame")
})
