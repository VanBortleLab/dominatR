test_that("plot_triangle works correctly", {

  # Sample data for testing
  data <- matrix(c(1:18), ncol=3, byrow=TRUE) |> as.data.frame()
  colnames(data) <- c("a", "b", "c")

  # Output is a dataframe
  result <- plot_triangle(data, output_table=TRUE)
  expect_true(is.data.frame(result))

  # Entropy and magnitude are normalized
  result <- plot_triangle(data, output_table=TRUE)
  expect_true(all(result$entropy >= 0 & result$entropy <= 1))
  expect_true(all(result$magnitude >= 0 & result$magnitude <= 1))
})
