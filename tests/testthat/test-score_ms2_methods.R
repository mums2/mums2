test_that("gnps_params works", {
  p <- gnps_params(frag_tolerance = 0.5)

  expect_equal(length(p), 2)
  expect_equal(p$tolerance, 0.5)

  expect_error(gnps_params())
})
