test_that("TPM normalization works correctly", {
  # Create a test data frame (gene expression matrix)
  test_data <- data.frame(
    sample1 = c(100, 200, 300, 400, 500),
    sample2 = c(50, 150, 250, 350, 450),
    sample3 = c(20, 120, 220, 320, 420)
  )

  # Create a gene length vector
  gene_length <- c(1000, 2000, 1500, 3000, 2500)

  # Calculate manual TPM
  df1 <- 1000 * test_data / gene_length
  expected_tpm <- as.matrix(sweep(df1, 2, colSums(df1) / 1e6, "/"))

  # Run TPM Normalization (without log2 transformation)
  result <- tpm_normalization(test_data, gene_length, log_trans = FALSE)

  # Expected to return a matrix
  expect_true(is.matrix(result))

  # Check that the number of rows and columns is maintained
  expect_equal(dim(result), dim(test_data))

  # Check that the TPM calculation is correct
  expect_equal(result, expected_tpm, tolerance = 1e-6)

  # Run the version with log2 transformation
  result_log <- tpm_normalization(test_data, gene_length, log_trans = TRUE)

  # Calculate the expected log2 transformation
  expected_log_tpm <- log2(expected_tpm + 1)

  # Check that the log2 transformation is correct
  expect_equal(result_log, expected_log_tpm, tolerance = 1e-6)
})
