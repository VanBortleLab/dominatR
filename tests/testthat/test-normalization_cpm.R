test_that("CPM normalization works correctly", {
  mat <- matrix(c(10,20,30,  40,50,60), nrow=2, byrow=FALSE)

  mat_cpm <- cpm_normalization(mat, log_trans=FALSE)

  # Check type
  expect_true(is.matrix(mat_cpm))
  expect_equal(dim(mat_cpm), dim(mat))

  col_sums <- colSums(mat, na.rm=TRUE)
  col_millions <- col_sums / 1e6
  expected <- sweep(mat, 2, col_millions, FUN="/")

  expect_equal(mat_cpm, expected, tolerance=1e-6)
})


test_that("cpm_normalization on matrix returns CPM matrix with log2 transform", {
  mat <- matrix(c(10,20,30,  40,50,60), nrow=2, byrow=TRUE)
  mat_cpm_log <- cpm_normalization(mat, log_trans=TRUE)

  # 1) first compute CPM
  col_sums <- colSums(mat)
  col_millions <- col_sums / 1e6
  cpm_vals <- sweep(mat, 2, col_millions, FUN="/")
  # 2) then log2(x+1)
  expected <- log2(cpm_vals + 1)

  expect_true(is.matrix(mat_cpm_log))
  expect_equal(dim(mat_cpm_log), dim(mat))
  expect_equal(mat_cpm_log, expected, tolerance=1e-6)
})

test_that("cpm_normalization fails for non-numeric input", {
  df_nonnum <- data.frame(x=c("a","b"), y=c("c","d"))
  expect_error(
    cpm_normalization(df_nonnum),
    "Input must be numeric"
  )
})

test_that("cpm_normalization on SummarizedExperiment overwrites default assay if new_assay_name not provided", {
  # Make a simple matrix
  mat <- matrix(c(10,20, 30,40), nrow=2, byrow=TRUE)
  colnames(mat) <- c("sample1", "sample2")
  rownames(mat) <- c("gene1", "gene2")
  se <- SummarizedExperiment(assays = list(counts=mat))

  se_cpm <- cpm_normalization(se)

  # Should overwrite "counts" with CPM
  cpm_mat <- assay(se_cpm, "counts")
  expect_true(is.matrix(cpm_mat))

  # Check that it's correct
  col_sums <- colSums(mat)
  col_millions <- col_sums / 1e6
  expected <- sweep(mat, 2, col_millions, FUN="/")
  expect_equal(cpm_mat, expected, tolerance=1e-6)
})

test_that("cpm_normalization on SummarizedExperiment creates new assay if new_assay_name is provided", {
  mat <- matrix(c(10,20,30, 40,50,60), nrow=2, byrow=TRUE)
  se <- SummarizedExperiment(assays = list(counts=mat))

  se_cpm <- cpm_normalization(se, new_assay_name="cpm_assay")

  # Original "counts" should remain unchanged
  expect_equal(assay(se_cpm, "counts"), mat)

  # The new assay "cpm_assay" should hold the normalized matrix
  expect_true("cpm_assay" %in% assayNames(se_cpm))

  cpm_mat <- assay(se_cpm, "cpm_assay")
  col_sums <- colSums(mat)
  col_millions <- col_sums / 1e6
  expected <- sweep(mat, 2, col_millions, FUN="/")
  expect_equal(cpm_mat, expected, tolerance=1e-6)
})

test_that("cpm_normalization handles log_trans=TRUE on SummarizedExperiment", {
  mat <- matrix(c(1,10, 100, 1000), nrow=2, byrow=TRUE)
  se <- SummarizedExperiment(assays = list(counts=mat))

  se_log <- cpm_normalization(se, log_trans=TRUE)
  # Overwrites the "counts" assay with log2-CPM
  log_mat <- assay(se_log, "counts")

  col_sums <- colSums(mat)
  col_millions <- col_sums / 1e6
  cpm_val <- sweep(mat, 2, col_millions, FUN="/")
  expected <- log2(cpm_val + 1)
  expect_equal(log_mat, expected, tolerance=1e-6)
})
