test_that("Assess MTLR_pred_model_f", {
  data(Train_Clin, package = "CPSM")
  data(Test_Clin, package = "CPSM")
  data(Key_Clin_feature_list, package = "CPSM")
  Result_Model_Type1 <- MTLR_pred_model_f(train_clin_data = Train_Clin,
                                          test_clin_data = Test_Clin,
                                          Model_type = 1,
                                          train_features_data = Train_Clin,
                                          test_features_data = Test_Clin,
                                          Clin_Feature_List = Key_Clin_feature_list,
                                          surv_time = "OS_month",
                                          surv_event = "OS",
                                          nfolds = 5)
    mean_median_survival_time_data <- Result_Model_Type1$mean_median_survival_time_data
    expect_s3_class(mean_median_survival_time_data, "data.frame")
})
