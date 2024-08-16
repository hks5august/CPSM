test_that("Assess_data_process", {
  data(Example_TCGA_LGG_FPKM_data, package = "SPM")
  New_data <- data_process_f(Example_TCGA_LGG_FPKM_data, col_num=20, surv_time = "OS.time")
  expect_s3_class(New_data, "data.frame")
  })
