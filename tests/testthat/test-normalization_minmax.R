test_that("minmax_normalization on matrix (one column) scales correctly", {
  mat <- matrix(c(10, 20, 30), nrow=3)
  # Single column => min=10, max=30 => range=20
  # new_min=0, new_max=1 => y=(x-10)/20
  scaled <- minmax_normalization(mat, new_min=0, new_max=1)
  expect_true(is.matrix(scaled))
  expect_equal(dim(scaled), dim(mat))

  # row1 => (10-10)/20=0, row2 =>(20-10)/20=0.5, row3 =>(30-10)/20=1
  expect_equal(scaled[,1], c(0, 0.5, 1))
})

test_that("minmax_normalization on matrix (multiple columns) with default range [0..1]", {
  mat <- matrix(c(1,2,3,  4,5,6), nrow=2, byrow=TRUE)

  scaled <- minmax_normalization(mat)
  expect_true(is.matrix(scaled))
  expect_equal(dim(scaled), c(2,3))

  # check a couple values
  # col1 row1 => (1-1)/3=0, row2 => (4-1)/3=1
  expect_equal(scaled[1,1], 0)
  expect_equal(scaled[2,1], 1)

  # col2 row1 => (2-2)/3=0, row2 => (5-2)/3=1
  expect_equal(scaled[1,2], 0)
  expect_equal(scaled[2,2], 1)

  # col3 row1 => (3-3)/3=0, row2 => (6-3)/3=1
  expect_equal(scaled[1,3], 0)
  expect_equal(scaled[2,3], 1)
})


test_that("minmax_normalization can scale to a user-specified range, e.g. [2..3]", {
  mat <- matrix(c(10,20,30), nrow=3)
  scaled <- minmax_normalization(mat, new_min=2, new_max=3)

  expect_equal(as.vector(scaled), c(2, 2.5, 3))
})

test_that("minmax_normalization fails on non-numeric input", {
  df <- data.frame(x=c("a","b"), y=c("c","d"))
  expect_error(minmax_normalization(df),
               "data is not numeric",
               fixed=TRUE)
})

test_that("minmax_normalization on SummarizedExperiment overwrites existing assay by default", {
  mat <- matrix(c(1,2, 3,4), nrow=2)
  se <- SummarizedExperiment(assays=list(counts=mat))
  se_scaled <- minmax_normalization(se)
  # overwrote "counts"
  expect_true(is(se_scaled, "SummarizedExperiment"))
  expect_equal(assayNames(se_scaled), "counts")

  # check scaled data
  scaled_assay <- assay(se_scaled, "counts")
  expect_true(is.matrix(scaled_assay))

  # check dimension
  expect_equal(dim(scaled_assay), dim(mat))

  expect_equal(scaled_assay[1,], c(0,0))
  expect_equal(scaled_assay[2,], c(1,1))
})

test_that("minmax_normalization on SummarizedExperiment can store in a new assay if new_assay_name is provided", {
  mat <- matrix(c(10,20,30,  40,50,60), nrow=2, byrow=TRUE)
  se <- SummarizedExperiment(assays=list(counts=mat))
  se2 <- minmax_normalization(se, new_assay_name="scaled_counts", new_min=-1, new_max=1)

  # original 'counts' should remain unchanged
  expect_true("counts" %in% assayNames(se2))
  expect_equal(assay(se2, "counts"), mat)

  # new assay
  expect_true("scaled_counts" %in% assayNames(se2))
  scaled_assay <- assay(se2, "scaled_counts")
  # check dimension
  expect_equal(dim(scaled_assay), c(2,3))

  expect_equal(scaled_assay[1,1], -1,  tolerance=1e-6)
  expect_equal(scaled_assay[2,1],  1,  tolerance=1e-6)
})
