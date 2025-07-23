test_that("dist_ms2 works with gnps parameters", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  dist <- dist_ms2(dat, 0.3, 2, gnps_params(0.5), min_peaks = 0)
  expect_s3_class(dist, "data.frame") 
  expect_equal(nrow(dist), 383)
})

test_that("dist_ms2 works with spectral entrophy parameters", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  dist <- dist_ms2(dat, 0.3, 2, spec_entropy_params(), min_peaks = 0)
  expect_s3_class(dist, "data.frame")
  expect_equal(nrow(dist), 357)

  dist <- dist_ms2(dat, 0.3, 2, spec_entropy_params(clean_spectra = FALSE),
  min_peaks = 0)
  expect_equal(nrow(dist), 591)

  dist <- dist_ms2(dat, 0.3, 2, spec_entropy_params(clean_spectra = FALSE, weighted = FALSE),
  min_peaks = 0)
  expect_equal(nrow(dist), 964)
})
