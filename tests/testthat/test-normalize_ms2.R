test_that("square root normalization works", {
  input <- c(20, 2, 103, 6)
  expected <- c(0.3907323, 0.1235604, 0.8867128, 0.2140129)

  out <- normalize_ms2(input, method = "sqrt")

  expect_equal(round(out, 7), expected)
})

test_that("scale normalization works", {
  input <- c(20, 2, 103, 6)
  expected <- c(0.1941748, 0.0194175, 1.0000000, 0.0582524)

  out <- normalize_ms2(input, method = "scale")

  expect_equal(round(out, 7), expected)
})
