test_that("plot_circle_frequency returns expected ggplot and data.frame structure", {
  # Create a fake 'circle' list whose second element has Entropy and Factor
  df_circle2 <- data.frame(
    Entropy = c(0.1, 0.2, 0.8, 0.9),
    Factor  = factor(c("X", "X", "Y", "Y"))
  )
  circle_obj <- list(NULL, df_circle2)

  result <- plot_circle_frequency(
    n           = 2,
    circle      = circle_obj,
    single      = TRUE,
    legend      = FALSE,
    numb_columns= 1,
    filter_class= NULL,
    point_size  = 1
  )

  expect_type(result, "list")
  expect_length(result, 2)
  expect_s3_class(result[[1]], "ggplot")
  expect_s3_class(result[[2]], "data.frame")

  required_cols <- c("bin", "Factor", "n", "proportion")
  expect_true(all(required_cols %in% names(result[[2]])))
})

test_that("plot_circle_frequency computes counts and proportions correctly", {
  # All Entropy values fall into a single bin when n = 1
  df_circle2 <- data.frame(
    Entropy = c(0.0, 0.0, 1.0, 1.0),
    Factor  = factor(c("A", "A", "B", "B"))
  )
  circle_obj <- list(NULL, df_circle2)

  result <- plot_circle_frequency(n = 2, circle = circle_obj, single = TRUE)
  df2 <- result[[2]]

  # Should have two row per Factor
  expect_equal(nrow(df2), 4)
  expect_setequal(as.character(df2$Factor), c("A", "B"))

  # Each factor has count = 2 and proportion = 1
  expect_true(all(df2$n %in% c(2,0)))
  expect_true(all(df2$proportion %in% c(1,0)))
})

test_that("plot_circle_frequency honors filter_class argument", {
  df_circle2 <- data.frame(
    Entropy = c(0.0, 0.0, 1.0, 1.0),
    Factor  = factor(c("A", "A", "B", "B"))
  )
  circle_obj <- list(NULL, df_circle2)

  result <- plot_circle_frequency(
    n            = 2,
    circle       = circle_obj,
    filter_class = "A"
  )
  df3 <- result[[2]]

  expect_equal(nrow(df3), 2)
  expect_equal(unique(df3$Factor), "A")
})
