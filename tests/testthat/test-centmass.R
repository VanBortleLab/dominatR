test_that("centmass on data.frame with default coords works as expected", {
  df <- data.frame(
    A = c(10, 0, 30),
    B = c( 5, 5,  0),
    C = c( 0, 0, 10)
  )

  df_res <- centmass(df)

  expect_true(is.data.frame(df_res))
  expect_equal(names(df_res), c("comx", "comy"))

  expect_equal(df_res$comx[1], 1/3,      tolerance=1e-6)
  expect_equal(df_res$comy[1], 0,        tolerance=1e-6)
  expect_equal(df_res$comx[2], 1,        tolerance=1e-6)
  expect_equal(df_res$comy[2], 0,        tolerance=1e-6)
  expect_equal(df_res$comx[3], 0.125,    tolerance=1e-6)

  # sqrt(3)/8 is about 0.21650635
  expect_equal(df_res$comy[3], sqrt(3)/8, tolerance=1e-7)
})

test_that("centmass on data.frame with custom x_coord, y_coord", {
  # Suppose we have 4 columns
  df <- data.frame(
    W = c(1, 2, 3),
    X = c(2, 4, 6),
    Y = c(3, 6, 9),
    Z = c(4, 8, 12)
  )

  # let's pick arbitrary coords
  x_coord <- c(1, 10, 100, 1000)
  y_coord <- c(2, 20, 200, 2000)
  df_res <- centmass(df, x_coord = x_coord, y_coord = y_coord)


  # For row1 sum=10
  expect_equal(df_res$comx[1], 432.1, tolerance=1e-6)
  expect_equal(df_res$comy[1], 864.2, tolerance=1e-6)

  # For row2 sum=20 => W=2, X=4, Y=6, Z=8
  expect_equal(df_res$comx[2], 432.1, tolerance=1e-6)
  expect_equal(df_res$comy[2], 864.2, tolerance=1e-6)

  # row3 sum=30 => W=3, X=6, Y=9, Z=12
  expect_equal(df_res$comx[3], 432.1, tolerance=1e-6)
  expect_equal(df_res$comy[3], 864.2, tolerance=1e-6)
})


test_that("centmass on data.frame fails if coord length mismatches ncol(data)", {
  df <- data.frame(A=1:3, B=4:6)
  # df has 2 columns, but x_coord,y_coord have length=3
  expect_error(centmass(df,
                        x_coord = c(0,1,0.5),
                        y_coord = c(0,0,sqrt(3)/2)),
               "Length of x_coord or y_coord does not match")
})

test_that("centmass on SummarizedExperiment works (default coords)", {
  mat <- matrix(c(10,0,30, 5,5,0, 0,0,10),
                nrow=3, byrow=FALSE)
  rownames(mat) <- paste0("gene", 1:3)
  colnames(mat) <- paste0("sample", 1:3)
  se <- SummarizedExperiment(assays = list(counts=mat))

  se_com <- centmass(se)

  # we expect rowData to have comx, comy
  expect_true(all(c("comx", "comy") %in% colnames(rowData(se_com))))
  expect_equal(unname(rowData(se_com)$comx[1]), 1/3, tolerance=1e-5)
  expect_equal(unname(rowData(se_com)$comy[1]), 0,   tolerance=1e-5)
  expect_equal(unname(rowData(se_com)$comx[2]), 1,   tolerance=1e-5)
  expect_equal(unname(rowData(se_com)$comy[2]), 0,   tolerance=1e-5)
  expect_equal(unname(rowData(se_com)$comx[3]), 0.125,      tolerance=1e-5)
  expect_equal(unname(rowData(se_com)$comy[3]), 0.216506 ,  tolerance=1e-5)
})

test_that("centmass on SummarizedExperiment fails if coordinate length mismatches ncol(assay)", {
  mat <- matrix(1:4, nrow=2)
  colnames(mat) <- c("s1","s2")
  se <- SummarizedExperiment(assays=list(counts=mat))
  # length mismatch: 2 col, but 3 coords
  expect_error(centmass(se, x_coord=c(0,1,2), y_coord=c(0,0,0)),
               "Length of x_coord or y_coord does not match")
})

test_that("centmass on SummarizedExperiment fails if no assay is present", {
  # SummarizedExperiment with no assay
  se_empty <- SummarizedExperiment()
  expect_error(centmass(se_empty),
               "No assays found",
               fixed=TRUE)

})
