test_that("dist_ms2_cpp works with gnps parameters", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  dist <- dist_ms2(dat, 0.3, 2, gnps_params(0.5), min_peaks = 0)
  expect_s3_class(dist, "data.frame")
  expect_equal(nrow(dist), 383)
})

test_that("dist_ms2_cpp works with spectral entrophy parameters", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  dist <- dist_ms2(dat, 0.3, 2, spec_entropy_params(), min_peaks = 0)
  expect_s3_class(dist, "data.frame")
  expect_equal(nrow(dist), 357)
})
