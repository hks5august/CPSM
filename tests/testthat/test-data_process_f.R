test_that("Assess_data_process", {
  data(Example_TCGA_LGG_FPKM_data, package = "CPSM")
  New_data <- data_process_f(assays(Example_TCGA_LGG_FPKM_data)$expression,
                             col_num=20, surv_time = "OS.time")
  expect_s3_class(New_data, "data.frame")
  })


