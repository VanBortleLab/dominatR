test_that("rpkm_normalization on matrix (no log) with user-supplied gene_length", {
  mat <- matrix(c(100, 200, 300,
                  400, 500, 600),
                nrow = 2, byrow = TRUE)
  rownames(mat) <- c("gene1", "gene2")
  colnames(mat) <- c("sample1", "sample2", "sample3")

  gene_len <- c(1000, 2000)  # both genes 1kb

  # Calculate RPKM
  res_rpkm <- rpkm_normalization(mat, gene_length = gene_len, log_trans = FALSE)

  # First gene, sample1:
  # count = 100, lib_size = 500, gene_length = 1000
  # RPKM = (100 * 1e9) / (500 * 1000) = 200000
  expect_true(is.matrix(res_rpkm))
  expect_equal(dim(res_rpkm), dim(mat))
  expect_equal(round(res_rpkm[1,1], 2), 200000.00, tolerance = 1e-2)
  # gene1, sample2: RPKM = (200 * 1e9) / (700 * 1000) = 285714.29
  expect_equal(round(res_rpkm[1,2], 2), 285714.29, tolerance = 1e-2)
  # gene2, sample2: RPKM = (500 * 1e9) / (700 * 2000) = 357142.86
  expect_equal(round(res_rpkm[2,2], 2), 357142.86, tolerance = 1e-2)
})


test_that("rpkm_normalization on SummarizedExperiment retrieves gene_length from rowData(x)$gene_length", {
  mat <- matrix(c(10, 20, 30, 40, 50, 60), nrow = 2, byrow = TRUE)
  rownames(mat) <- c("gene1", "gene2")
  colnames(mat) <- c("sampleA", "sampleB", "sampleC")
  se <- SummarizedExperiment(assays = list(counts = mat))
  rowData(se)$gene_length <- c(1000, 2000)

  se_rpkm <- rpkm_normalization(se, log_trans = FALSE)

  rpkm_mat <- assay(se_rpkm, "counts")
  expect_true(is.matrix(rpkm_mat))
  expect_equal(round(rpkm_mat[1,1], 2), 200000.00, tolerance = 1e-2)
})

test_that("rpkm_normalization on SummarizedExperiment can store in a new assay if new_assay_name is given", {
  mat <- matrix(1:6, nrow = 2)
  se <- SummarizedExperiment(assays = list(mycounts = mat))
  rowData(se)$gene_length <- c(100, 200)

  se_new <- rpkm_normalization(se, assay_name = "mycounts", new_assay_name = "RPKM", log_trans = TRUE)

  expect_equal(assay(se_new, "mycounts"), mat)
  expect_true("RPKM" %in% assayNames(se_new))
})

test_that("rpkm_normalization on SummarizedExperiment fails with incorrect or missing gene_length", {
  mat <- matrix(1:6, nrow = 2)
  se <- SummarizedExperiment(assays = list(counts = mat))

  expect_error(rpkm_normalization(se), "No 'gene_length' column found in rowData\\(x\\)")

  rowData(se)$gene_length <- c("a", "b")
  expect_error(rpkm_normalization(se), "'gene_length' in rowData\\(x\\) must be numeric")
})
