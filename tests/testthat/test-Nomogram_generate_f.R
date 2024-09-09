test_that("Assess Nomogram_generate_f", {
  # Load necessary data
  data(Train_Data_Nomogram_input, package = "CPSM")
  data(feature_list_for_Nomogram, package = "CPSM")

  # Run the Nomogram_generate_f function
  Result_Nomogram <- Nomogram_generate_f(
    data = Train_Data_Nomogram_input,
    Feature_List = feature_list_for_Nomogram,
    surv_time = "OS_month",
    surv_event = "OS"
  )

  # Check that the result is a list
  expect_type(Result_Nomogram, "list")

  # Check that the list contains the correct named element
  expect_true("C_index_mat" %in% names(Result_Nomogram))

  # Verify that C_index_mat is a matrix with correct dimensions
  C_index_mat <- Result_Nomogram$C_index_mat
  expect_true(is.matrix(C_index_mat))
  expect_equal(dim(C_index_mat), c(1, 2))

  # Check that C_index_mat has the correct column names
  expect_equal(colnames(C_index_mat), c("Bias-corrected C-index", "C-index"))

  # Optionally, you can check the range of the c-index values
  expect_true(all(C_index_mat >= 0 & C_index_mat <= 1))
})
