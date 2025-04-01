
test_that("entropy on data.frame row-normalizes and computes Entropy", {
  # Prepare test data
  df <- data.frame(
    A = c(10, 0, 30),
    B = c(5,  5,  0),
    C = c(0,  0,  10)
  )

  df_e <- entropy(df)

  expected_norm <- data.frame(
    A = c(2/3, 0, 3/4),
    B = c(1/3, 1, 0),
    C = c(0,0, 1/4 ))




  expect_equal(df_e[, c("A","B","C")],
               expected_norm,
               tolerance = 1e-6)

  # Check that Entropy column exists
  expect_true("Entropy" %in% names(df_e))

  expect_equal(df_e$Entropy[1], 0.9183, tolerance=1e-3)
  expect_equal(df_e$Entropy[2], 0,      tolerance=1e-3)
  expect_equal(df_e$Entropy[3], 0.8113, tolerance=1e-3)
})

test_that("entropy on data.frame with no numeric columns returns unchanged or warns", {
  df_no_num <- data.frame(X=letters[1:3])
  # We expect a warning but return the same df
  expect_warning({
    df_res <- entropy(df_no_num)
    expect_equal(df_res, df_no_num)
  })
})




test_that("entropy on SummarizedExperiment normalizes assay and appends rowData$Entropy", {
  mat <- matrix(c(10,0,30, 5,5,0, 0,0,10), nrow=3, byrow=F)
  rownames(mat) <- paste0("gene", 1:3)
  colnames(mat) <- paste0("sample", 1:3)
  se <- SummarizedExperiment(assays=list(counts=mat))

  # Apply entropy
  se_e <- entropy(se, new_assay_name = 'Entropy')
  # The 'counts' assay should now be row-normalized
  norm_counts <- assay(se_e, "Entropy")

  expected_norm <- cbind(
      sample1 = c(2/3, 0, 3/4),
      sample2 = c(1/3, 1, 0),
      sample3 = c(0,0, 1/4 ))

  rownames(expected_norm) = paste0("gene", 1:3)

  expect_equal(norm_counts, expected_norm, tolerance=1e-6)

  # rowData should have Entropy
  expect_true("Entropy" %in% colnames(rowData(se_e)))
})

