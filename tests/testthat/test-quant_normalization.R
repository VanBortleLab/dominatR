test_that("Quantile normalization works correctly", {
  test_data <- data.frame(
    sample1 = c(10, 20, 30, 40, 50),
    sample2 = c(5, 15, 25, 35, 45),
    sample3 = c(2, 12, 22, 32, 42)
  )

  result <- quant_normalization(test_data, log_trans = FALSE)
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(test_data))
  result_log <- quant_normalization(test_data, log_trans = TRUE)
  expect_equal(result_log, log2(result + 1))
  expect_equal(rownames(result), rownames(test_data))
})
