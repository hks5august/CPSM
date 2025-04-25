test_that("multiplication works", {
data(Train_PI_data, package = "CPSM")
data(Test_PI_data, package = "CPSM")
data(Key_PI_list, package = "CPSM")

Results_Risk_group_Prediction<-  predict_survival_risk_group_f(selected_train_data = Train_PI_data,
                                                             selected_test_data = Test_PI_data,
                                                             Feature_List = Key_PI_list )  
  expect_s3_class(Results_Risk_group_Prediction$Test_results, "data.frame")
})
