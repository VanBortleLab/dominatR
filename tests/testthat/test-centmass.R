test_that("Center of mass works correctly", {

  ## Hypothethical dataframe
  data <- matrix(c(1, 1, 1), nrow = 1)

  ## Calculate values with the function
  result <- centmass(data)

  ### This is what we expect for each coordinate
  expected_x <- (0*1 + 1*1 + 0.5*1) / 3
  expected_y <- (0*1 + 0*1 + (sqrt(3)/2)*1) / 3

  ## Equal Values
  expect_equal(result$comx, expected_x)
  expect_equal(result$comy, expected_y)

  ### 0 Values
  data <- matrix(c(0, 0, 0), nrow = 1)

  ### Should be 0
  result <- centmass(data)
  # Expect NaN or some defined behavior, depends on how you want to handle this case
  expect_true(all(is.nan(result$comx)))
  expect_true(all(is.nan(result$comy)))

})
