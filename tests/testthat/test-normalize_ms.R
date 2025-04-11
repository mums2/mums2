test_that("relative_abundance works", {
  raw_values <- c(0.00, 1.00, 0.00, 5.00, 2.50, 1.50)

  norm_values <- relative_abundance(raw_values)
  expect_equal(norm_values, c(0, 0.1, 0, 0.5, 0.25, 0.15))
})


