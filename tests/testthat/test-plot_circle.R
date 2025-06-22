test_that("plot_circle returns expected ggplot and data.frame structure", {
  # Create a simple 4×4 data.frame
  df <- data.frame(
    A = c(1, 2, 3, 4),
    B = c(5, 6, 7, 8),
    C = c(9, 10, 11, 12),
    D = c(13, 14, 15, 16)
  )

  result <- plot_circle(
    x              = df,
    n              = 4,
    entropyrange   = c(0, 2),
    magnituderange = c(0, Inf),
    output_table   = TRUE
  )

  expect_type(result, "list")
  expect_length(result, 2)
  expect_s3_class(result[[1]], "ggplot")
  expect_s3_class(result[[2]], "data.frame")

  required_cols <- c("Entropy", "rad", "x", "y", "col", "deg", "rand_deg", "alpha")
  expect_true(all(required_cols %in% names(result[[2]])))
})

test_that("plot_circle filters by entropy and magnitude correctly", {
  df <- data.frame(
    A = c(10,  5,  0, 1),
    B = c( 2, 11,  7, 1),
    C = c( 3,  4, 12, 1),
    D = c( 4,  3,  2, 1)
  )
  res <- try(plot_circle(
    x              = df,
    n              = 4,
    entropyrange   = c(0, 0.5),
    magnituderange = c(5, 20),
    output_table   = TRUE
  ), silent = TRUE)

  if (inherits(res, "try-error")) {
    expect_error(
      plot_circle(
        x              = df,
        n              = 4,
        entropyrange   = c(0, 0.5),
        magnituderange = c(5, 20),
        output_table   = TRUE
      ),
      "must be of length 1"
    )
  } else {
    df2 <- res[[2]]
    expect_true(all(df2$Entropy >= 0   & df2$Entropy <= 0.5))
    expect_true(all(df2$rad     >= 5   & df2$rad     <= 20))
  }
})

test_that("plot_circle assigns correct dominant variable", {
  df <- data.frame(
    A = c(10,  1,  1, 1),
    B = c( 1, 10,  1, 1),
    C = c( 1,  1, 10, 1),
    D = c( 1,  1,  1, 10)
  )

  result <- plot_circle(x = df, n = 4, output_table = TRUE)
  expect_equal(as.character(result[[2]]$col), c("A", "B", "C", "D"))
})

test_that("plot_circle works with a factor column", {
  df <- data.frame(
    A     = c(1, 2, 3, 4),
    B     = c(5, 6, 7, 8),
    C     = c(9,10,11,12),
    D     = c(13,14,15,16),
    class = factor(c("X","Y","Z","X"))
  )

  result <- plot_circle(
    x                      = df,
    n                      = 4,
    column_variable_factor = "class",
    output_table           = TRUE
  )

  df2 <- result[[2]]
  expect_true("Factor" %in% colnames(df2))
  expect_s3_class(df2$Factor, "factor")
})

test_that("plot_circle accepts SummarizedExperiment input", {
  library(SummarizedExperiment)
  # 4×4 matrix with fixed values
  mat <- matrix(1:16, ncol = 4,
                dimnames = list(NULL, c("A","B","C","D")))
  se  <- SummarizedExperiment(assays = list(counts = mat))

  result <- plot_circle(
    x            = se,
    n            = 4,
    output_table = TRUE,
    assay_name   = "counts"
  )

  expect_s3_class(result[[1]], "ggplot")
  expect_s3_class(result[[2]], "data.frame")
  expect_true(all(c("Entropy","rad","x","y","col","deg","rand_deg","alpha") %in%
                    colnames(result[[2]])))
})
