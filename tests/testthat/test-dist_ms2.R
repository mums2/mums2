test_that("dist_ms2 works with gnps parameters", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  dist <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5), min_peaks = 0)
  expect_s3_class(dist, "data.frame")
  expect_equal(nrow(dist), 383)
})

test_that("dist_ms2 works with spectral entropy parameters", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  dist <- dist_ms2(dat, 0.3, 2, spec_entropy_params(), min_peaks = 0)
  expect_s3_class(dist, "data.frame")
  expect_equal(nrow(dist), 357)

  dist <- dist_ms2(dat, 0.3, 2, spec_entropy_params(clean_spectra = FALSE),
                   min_peaks = 0)
  expect_equal(nrow(dist), 591)

  dist <- dist_ms2(dat, 0.3, 2, spec_entropy_params(clean_spectra = FALSE,
                                                    weighted = FALSE),
                   min_peaks = 0)
  expect_equal(nrow(dist), 964)
})


test_that("dist_ms2 error conditions are working as expected", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  expect_error(dist_ms2("dat", 0.3, 2, spec_entropy_params(), min_peaks = 0),
              "The mass_data object must be created using")
  expect_error(dist_ms2(dat, "0.3", 2, spec_entropy_params(), min_peaks = 0),
              "cutoff must be numeric")
  expect_error(dist_ms2(dat, 0.3, "2", spec_entropy_params(), min_peaks = 0),
              "precursor_threshold must be a numeric")
  expect_error(dist_ms2(dat, 0.3, 2, c(), min_peaks = 0),
              "score_params must be created")
  expect_error(dist_ms2(dat, 0.3, 2, spec_entropy_params(), min_peaks = ""),
              "min_peaks must be a numeric")
  expect_error(dist_ms2(dat, 0.3, 2, spec_entropy_params(), number_of_threads = ""),
              "number_of_threads must be a numeric")

})

test_that("fails when there are no ms2 matches", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  dat$ms2_matches <- data.frame()
  expect_error(dist_ms2(dat, 0.3, 2, spec_entropy_params(), min_peaks = 0),
               "Cannot calculate distances")
})
