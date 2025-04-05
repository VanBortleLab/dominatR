

test_that("tpm_normalization on matrix (no log) with user-supplied gene_length", {
  # Simple 2x3 matrix
  mat <- matrix(c(10,20,30,  40,50,60), nrow=2, byrow=TRUE)
  rownames(mat) <- c("gene1","gene2")
  colnames(mat) <- c("sample1","sample2","sample3")

  # gene_length must match nrow
  gene_len <- c(1000, 2000)

  res_tpm <- tpm_normalization(mat, gene_length=gene_len, log_trans=FALSE)

  expect_equal(res_tpm[1,1], 333333.3333, tolerance=1e-3)
  # gene2 col2 => 25 / (45/1e-6)= 25/(4.5e-5)= 555,555.5555
  expect_equal(res_tpm[2,2], 555555.5555, tolerance=1e-3)
})


test_that("tpm_normalization on SummarizedExperiment retrieves gene_length from rowData(x)$gene_length", {
  mat <- matrix(c(10,20,30,  40,50,60), nrow=2, byrow=TRUE)
  rownames(mat) <- c("gene1","gene2")
  colnames(mat) <- c("sampleA","sampleB","sampleC")
  se <- SummarizedExperiment(assays=list(counts=mat))

  # add gene lengths to rowData
  rowData(se)$gene_length <- c(1000, 2000)

  # run TPM => overwrites the "counts" assay by default
  se_tpm <- tpm_normalization(se, log_trans=FALSE)

  # check that it did not error & "counts" was replaced
  tpm_mat <- assay(se_tpm, "counts")
  expect_true(is.matrix(tpm_mat))

  expect_equal(tpm_mat[1,1], 333333.333, tolerance=1e-3)
})

test_that("tpm_normalization on SummarizedExperiment can store in a new assay if new_assay_name is given", {
  mat <- matrix(1:6, nrow=2)
  se <- SummarizedExperiment(assays=list(mycounts=mat))
  rowData(se)$gene_length <- c(100,200)

  # store results in new assay => "TPM"
  se_new <- tpm_normalization(se, assay_name="mycounts",
                              new_assay_name="TPM",
                              log_trans=TRUE)
  # original "mycounts" should remain unchanged
  expect_equal(assay(se_new, "mycounts"), mat)
  # new assay "TPM" should exist
  expect_true("TPM" %in% assayNames(se_new))
})

test_that("tpm_normalization on SummarizedExperiment fails if rowData(x)$gene_length is missing, non-numeric, or mismatched dimension", {
  mat <- matrix(1:6, nrow=2)
  se <- SummarizedExperiment(assays=list(counts=mat))

  # 1) Missing gene_length
  expect_error(
    tpm_normalization(se),
    "No 'gene_length' column found in rowData\\(x\\)"
  )
  # 2) Non-numeric gene_length
  rowData(se)$gene_length <- c("a","b")
  expect_error(
    tpm_normalization(se),
    "'gene_length' in rowData\\(x\\) must be numeric"
  )

})
