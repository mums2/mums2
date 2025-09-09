test_that("modified_cosine_params works", {
  p <- modified_cosine_params(frag_tolerance = 0.5)

  expect_equal(length(p), 2)
  expect_equal(p$tolerance, 0.5)

  expect_error(modified_cosine_params())
})

test_that("spec_entropy_params works", {
  p <- spec_entropy_params()

  expect_equal(length(p), 9)
  expect_equal(p$clean_spectra, TRUE)
})
