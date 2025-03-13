test_that("minmax normalization works correctly", {

   ### An Expected Output, values within a range
    x <- c(1, 2, 3, 4, 5)
    result <- minmax_normalization(x, 0, 1)
    expect_true(all(result >= 0 & result <= 1))

  ### Standard min max normalization for 3 known values
    x <- c(10, 20, 30)
    known_output <- c(0, 0.5, 1)
    result <- minmax_normalization(x, 0, 1)
    expect_equal(result, known_output)
  })
