test_that("multiplication works", {
  data(Train_Norm_data, package = "CPSM")
  data(Test_Norm_data, package = "CPSM")
  Result_PI <- Lasso_PI_scores_f(train_data = Train_Norm_data,
                                 test_data = Test_Norm_data,
                                 nfolds=5,
                                 n_repeats = 10,    
    				 alpha = 1,            
				 col_num=21,
                                 freq_threshold = 0.6,
                                 surv_time = "OS_month",
                                 surv_event = "OS")

# Check that key outputs exist and are data.frames
  expect_s3_class(Result_PI$Train_Lasso_key_variables, "data.frame")
  expect_s3_class(Result_PI$Train_PI_data, "data.frame")
  expect_s3_class(Result_PI$Test_PI_data, "data.frame")
  expect_s3_class(Result_PI$Stability_table, "data.frame")

  # Optionally: check Stability_table has reasonable values
  expect_true(all(Result_PI$Stability_table$Selection_Frequency >= 0 &
                  Result_PI$Stability_table$Selection_Frequency <= 1))
})
