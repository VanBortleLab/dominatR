test_that("quantile_normalization works on matrix (no log)", {
  mat <- matrix(c(1, 4, 3, 2, 5, 6), nrow = 3, ncol = 2)
  rownames(mat) <- paste0("gene", 1:3)
  colnames(mat) <- paste0("sample", 1:2)

  norm_mat <- quantile_normalization(mat, log_trans = FALSE)

  expect_true(is.matrix(norm_mat))
  expect_equal(dim(norm_mat), dim(mat))
  expect_equal(round(norm_mat[1,1], 2), round(norm_mat[1,2], 2), tolerance = 1e-6)  # same quantile
})

test_that("quantile_normalization works on SummarizedExperiment", {
  mat <- matrix(c(1, 4, 3, 2, 5, 6), nrow = 3, ncol = 2)
  rownames(mat) <- paste0("gene", 1:3)
  colnames(mat) <- paste0("sample", 1:2)
  se <- SummarizedExperiment(assays = list(counts = mat))

  se_qn <- quantile_normalization(se)
  norm_mat <- assay(se_qn, "counts")

  expect_true(is.matrix(norm_mat))
  expect_equal(dim(norm_mat), dim(mat))
})

test_that("quantile_normalization stores in new assay if new_assay_name is provided", {
  mat <- matrix(1:6, nrow = 3)
  se <- SummarizedExperiment(assays = list(mycounts = mat))

  se_new <- quantile_normalization(se, assay_name = "mycounts", new_assay_name = "QNorm")

  expect_equal(assay(se_new, "mycounts"), mat)
  expect_true("QNorm" %in% assayNames(se_new))
})

test_that("quantile_normalization fails on non-numeric data.frame or matrix", {
  df <- data.frame(a = c("a", "b", "c"), b = c("x", "y", "z"))
  expect_error(quantile_normalization(df), "Input must be numeric")

  mat <- matrix(letters[1:6], nrow = 2)
  expect_error(quantile_normalization(mat), "Input must be numeric")
})
