test_that("Assess_data_process", {
  data(Example_TCGA_LGG_FPKM_data, package = "CPSM")
  combined_df <- cbind(as.data.frame(colData(Example_TCGA_LGG_FPKM_data))
                       [, -ncol(colData(Example_TCGA_LGG_FPKM_data))],
                       t(as.data.frame(assay(Example_TCGA_LGG_FPKM_data,
                                             "expression"))))
  New_data <- data_process_f(combined_df, col_num=20, surv_time = "OS.time")
  expect_s3_class(New_data, "data.frame")
  })


