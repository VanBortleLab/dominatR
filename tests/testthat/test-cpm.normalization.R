test_that("CPM normalization works correctly", {
  test_data <- data.frame(
    sample1 = c(100, 200, 300, 400, 500),
    sample2 = c(50, 150, 250, 350, 450),
    sample3 = c(20, 120, 220, 320, 420)
  )
  expected_cpm <- sweep(test_data, 2, colSums(test_data) / 1e6, "/")
  result <- cpm.normalization(test_data, log_trans = FALSE)
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(test_data))
  expect_equal(result, expected_cpm, tolerance = 1e-6)
  result_log <- cpm.normalization(test_data, log_trans = TRUE)
  expected_log_cpm <- log2(expected_cpm + 1)
  expect_equal(result_log, expected_log_cpm, tolerance = 1e-6)
})
