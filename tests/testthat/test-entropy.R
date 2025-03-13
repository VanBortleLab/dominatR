test_that("Entropy function works", {
  data <- tibble::tibble(a = c(1,2,3), b = c(4,5,6))
  result <- entropy(data)
  expect_true("Entropy" %in% colnames(result))

  ### Max Entropy is equal to log2(N) if all variables contribute equally
  data <- tibble::tibble(a = c(1,1,1), b = c(1,1,1))
  result <- entropy(data)


  values = c(1,1,1)

  expect_equal(result$Entropy, values)

})
