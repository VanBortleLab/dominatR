
test_that("plot_rope returns an augmented data frame for ordinary data", {

  ## two simple numeric columns
  tbl <- data.frame(A = c(10,  5,  0),
                    B = c( 2, 11,  7))

  ## call the function with output_table = TRUE
  res <- plot_rope(tbl,
                   push_text     = 1,
                   rope_width    = 1,
                   output_table  = TRUE)

  ##  --- structural checks ---------------------------------
  expect_s3_class(res, "data.frame")


  ## The function returns 7 columns inf output_table = TRUEF
  expect_equal(ncol(res), 2 + 5)
  expect_setequal(colnames(res),
                  c("a","b", "comx","comy",
                    "color","maxvalue","entropy"))

  red_dim <- length(res$color[ res$a > res$b ])
  blue_dim <- length(res$color[ res$a <= res$b ])

  expect_equal(res$color[ res$a > res$b ], rep("red", red_dim))
  expect_equal(res$color[ res$a <= res$b ], rep("blue", blue_dim))

})


test_that("plot_rope handles efficiently data according to an entropy range", {

  ## two simple numeric columns
  tbl <- data.frame(A = c(10,  5,  0),
                    B = c( 2, 11,  7))

  ## call the function with output_table = TRUE
  res <- plot_rope(tbl,
                   push_text     = 1,
                   rope_width    = 1,
                   entropyrange = c(0,0.2),
                   output_table  = TRUE)

  ##  --- structural checks ---------------------------------
  expect_s3_class(res, "data.frame")

  ## the function filters out values outside of the entropy range
  red_dim <- length(res$color[ res$entropy  < 0.2])
  blue_dim <- length(res$color[ res$entropy > 0.2])

  expect_equal(res$color[ res$entropy < 0.2 ], rep("blue", red_dim))
  expect_equal(res$color[ res$a > 0.2 ], rep("whitesmoke", blue_dim))

})



test_that("plot_rope handles efficiently data according to an entropy range and maxvalue range", {

  ## two simple numeric columns
  tbl <- data.frame(A = c(10,  5,  0, 100),
                    B = c( 2, 11,  7, 1000))

  ## call the function with output_table = TRUE
  res <- plot_rope(tbl,
                   push_text     = 1,
                   rope_width    = 1,
                   entropyrange = c(0,0.2),
                   maxvaluerange = c(0,50),
                   output_table  = TRUE)

  ##  --- structural checks ---------------------------------
  expect_s3_class(res, "data.frame")

  ## the function filters out values outside of the entropy range
  red_dim <- length(res$color[ res$entropy  < 0.2 & res$maxvalue < 100])
  blue_dim <- length(res$color[ res$entropy > 0.2 & res$maxvalue > 100])

  expect_equal(res$color[ res$entropy < 0.2 & res$maxvalue < 100],
               rep("blue", red_dim))
  expect_equal(res$color[ res$a > 0.2 & res$maxvalue > 100],
               rep("whitesmoke", blue_dim))

})


test_that("plot_rope works with SummarizedExperiment + column_name", {

  mat <- matrix(c(4,1,  6,9,   # two samples, three genes
                  3,8),
                nrow = 3, ncol = 2, byrow = FALSE,
                dimnames = list(paste0("g",1:3),
                                c("condA","condB")))
  se  <- SummarizedExperiment(assays = list(counts = mat))

  tbl <- plot_rope(se,
                   column_name   = c("condA","condB"),
                   assay_name    = "counts",
                   output_table  = TRUE)

  expect_s3_class(tbl, "data.frame")
  expect_equal(nrow(tbl), 3)
  expect_equal(tbl$a, unname(mat[,1]))
  expect_equal(tbl$b, unname(mat[,2]))

})


