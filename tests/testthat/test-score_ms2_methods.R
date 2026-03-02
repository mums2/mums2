test_that("modified_cosine_params works", {
  parameters <- modified_cosine_params(frag_tolerance = 0.5)

  expect_equal(length(parameters), 2)
  expect_equal(parameters$tolerance, 0.5)
  expect_s3_class(parameters, "parameters")
  expect_error(modified_cosine_params())
})

test_that("spec_entropy_params works", {
  parameters <- spec_entropy_params()

  expect_equal(length(parameters), 9)
  expect_s3_class(parameters, "parameters")
  expect_equal(parameters$clean_spectra, TRUE)
})
