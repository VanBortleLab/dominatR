
test_that("plot_triangle returns augmented data for matrix input", {

  mat <- matrix(c(10,  2,  3,
                  4,  15,  1,
                  7,   6, 12),
                ncol = 3, byrow = TRUE,
                dimnames = list(NULL, c("A","B","C")))

  res <- plot_triangle(mat,
                       entropyrange   = c(0, Inf),
                       maxvaluerange  = c(0, Inf),
                       output_table   = TRUE,
                       plotAll        = FALSE)

  ## structure
  expect_s3_class(res, "data.frame")
  expect_equal(ncol(res), 3 + 5)                # a,b,c + 5 extras
  expect_setequal(colnames(res),
                  c("a","b","c",
                    "Entropy","comx","comy",
                    "max_counts","color"))

  ## colour assignment
  expect_equal(res$color[1], "darkred")   # a is max
  expect_equal(res$color[2], "darkgreen") # b is max
  expect_equal(res$color[3], "darkblue")  # c is max
})


test_that("entropy / max filters turn colour to background", {

  mat <- matrix(c(100,1,1,
                  1,1,1), ncol = 3, byrow = TRUE)

  ## very tight max filter to exclude second row
  res <- plot_triangle(mat,
                       maxvaluerange = c(50, 200),
                       output_table  = TRUE)

  expect_equal(res$color[2], "whitesmoke")
  expect_true(res$color[1] != "whitesmoke")
})

test_that("plot_triangle handles SummarizedExperiment and column_name", {

  mat <- matrix(c(5, 9, 2,
                  2, 3, 8,
                  6, 2, 1),
                ncol = 3, byrow = TRUE,
                dimnames = list(paste0("g",1:3),
                                c("cond1","cond2","cond3")))
  se <- SummarizedExperiment(assays = list(cnt = mat))

  res <- plot_triangle(se,
                       column_name = c("cond1","cond2","cond3"),
                       assay_name  = "cnt",
                       output_table = TRUE)

  ## first row: cond3=2, cond1=5, cond2=9  => max=cond2 (darkgreen)
  expect_equal(res$color[1], "darkgreen")
})

test_that("plot_triangle errors on wrong column count or non-numeric data", {

  mat_bad <- matrix(1:4, ncol = 2)      # only two columns
  expect_error(plot_triangle(mat_bad, output_table = TRUE))

  df_char <- data.frame(A = letters[1:3],
                        B = letters[1:3],
                        C = letters[1:3])
  expect_error(plot_triangle(df_char, output_table = TRUE))
})
