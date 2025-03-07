test_that("Entropy function works", {
  data <- tibble::tibble(a = c(1,2,3), b = c(4,5,6))
  result <- entropy(data)
  expect_true("Entropy" %in% colnames(result))
})
