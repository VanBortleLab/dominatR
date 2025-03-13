test_that("RPKM normalization works correctly", {
  # Create a test data frame (gene expression matrix)
  test_data <- data.frame(
    sample1 = c(100, 200, 300, 400, 500),
    sample2 = c(50, 150, 250, 350, 450),
    sample3 = c(20, 120, 220, 320, 420)
  )

  # Create a gene length vector
  gene_length <- c(1000, 2000, 1500, 3000, 2500)

  # Calculate manual RPKM
  counts_per_million <- sweep(test_data, 2, colSums(test_data) / 1e6, "/")
  expected_rpkm <- sweep(counts_per_million * 1000, 1, gene_length, "/")

  # Run RPKM normalization (without log2 transformation)
  result <- rpkm.normalization(test_data, gene_length, log_trans = FALSE)

  # Expected to return a matrix(according to return)
  expect_true(is.matrix(result))

  # Check that the number of rows and columns is maintained
  expect_equal(dim(result), dim(test_data))

  # Check that RPKM calculation is correct
  expect_equal(result, expected_rpkm, tolerance = 1e-6)

  # Run the version with log2 transformation
  result_log <- rpkm.normalization(test_data, gene_length, log_trans = TRUE)

  # Calculate the expected log2 transformation
  expected_log_rpkm <- log2(expected_rpkm + 1)

  # Check that the log2 transformation is correct
  expect_equal(result_log, expected_log_rpkm, tolerance = 1e-6)
})
