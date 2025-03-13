test_that("Qentropy works correctly", {
  data <- tibble::tibble(class1 = c(4,2,3), class2 = c(1,1,9), class3 = c(0,0,0))
  result = entropy(data)

  result = Qentropy(result)

  ### For columns with 0 after Entropy calculation values should be infinite
  expect_equal(result$class3, rep(Inf, 3))

  ### Final dataframe should not have Entropy values
  expect_false('Entropy' %in% colnames(result))

})
