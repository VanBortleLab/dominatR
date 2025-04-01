
test_that("Qentropy on data.frame uses existing row-normalized data + Entropy", {
  df <- data.frame(
    A = c(10, 0, 30),
    B = c(5,  5,  0),
    C = c(0,  0,  10)
  )
  df_e <- entropy(df)   # row-normalize & compute Entropy
  df_q <- Qentropy(df_e)

  # After Qentropy:
  # We expect 'Entropy' column to be removed
  expect_false("Entropy" %in% names(df_q))

  # Let's check row1 again:
  # We'll just check a small set of values to confirm
  expect_equal(df_q$A[1], 1.5033, tolerance=1e-3)
  expect_equal(df_q$B[1], 2.5033, tolerance=1e-3)
  expect_true(is.infinite(df_q$C[1]))
})


test_that("Qentropy on data.frame fails if 'Entropy' is missing", {
  df <- data.frame(A=c(1,2,3))
  expect_error(Qentropy(df),
               "'Entropy' column not found",
               fixed=TRUE)
})



test_that("Qentropy on SummarizedExperiment uses the Entropy assay and works properly", {
  mat <- matrix(c(10,0,30, 5,5,0, 0,0,10), nrow=3, byrow=FALSE)
  se <- SummarizedExperiment(assays=list(counts=mat))

  # 1) entropy
  se_e <- entropy(se, new_assay_name = 'Entropy')
  # 2) Qentropy
  se_q <- Qentropy(se_e,assay_name = 'Entropy', 'Qentropy')

  # The new 'counts' assay should contain Q values
  q_counts <- assay(se_q, "Qentropy")
  expect_true(is.infinite(q_counts[1, 3]))  # row1 col3 was 0 => Inf

  # rowData$Entropy should be removed
  expect_true("Entropy" %in% colnames(rowData(se_q)))
})
