test_that("CPM normalization works correctly", {
  test_data <- data.frame(
    sample1 = c(100, 200, 300, 400, 500),
    sample2 = c(50, 150, 250, 350, 450),
    sample3 = c(20, 120, 220, 320, 420)
  )
  ## Our primary calculation - returns a matrix
  expected_cpm <- as.matrix(sweep(test_data, 2, colSums(test_data) / 1e6, "/"))

  ### Using our function
  result <- cpm_normalization(test_data, log_trans = FALSE)

  ## It should be a matrix
   expect_true(is.matrix(result))

  ## Same number of dimensions
  expect_equal(dim(result), dim(test_data))

  ## values should be equal
  expect_equal(result, expected_cpm, tolerance = 1e-6)

  ## Using the log feature.
  result_log <- cpm_normalization(test_data, log_trans = TRUE)

  ## Adding a count to our primary calculation
  expected_log_cpm <- log2(expected_cpm + 1)

  ## Should be equal
  expect_equal(result_log, expected_log_cpm, tolerance = 1e-6)
})
