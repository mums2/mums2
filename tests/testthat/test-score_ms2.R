test_that("score_ms2 works on class peaks_data", {

  dir <- "exttestdata"
  file1 <- "MSV000086713_32.csv"
  file2 <- "PSUMSMLS_Adenine.csv"

  scan_32 <- create_peaks_data(read.csv(test_path(dir, file1)),
                               136.0615)

  adenine <- create_peaks_data(read.csv(test_path(dir, file2)),
                               136.0620)

  score <- score_ms2(scan_32, adenine, gnps_params(frag_tolerance = 0.5))
  expect_equal(score, 0.81162792)
})

test_that("score_ms2 works on class peaks_data, method entropy", {

  dir <- "exttestdata"
  file1 <- "MSV000086713_32.csv"
  file2 <- "PSUMSMLS_Adenine.csv"

  scan_32 <- create_peaks_data(read.csv(test_path(dir, file1)),
                               136.0615)

  adenine <- create_peaks_data(read.csv(test_path(dir, file2)),
                               136.0620)

  e <- msentropy::msentropy_similarity(as.matrix(scan_32$spectra),
                                       as.matrix(adenine$spectra),
                                       ms2_tolerance_in_da = 0.02,
                                       ms2_tolerance_in_ppm = -1,
                                       clean_spectra = TRUE,
                                       min_mz = 0,
                                       max_mz = 1000,
                                       noise_threshold = 0.01,
                                       max_peak_num = 100,
                                       weighted = TRUE)

  score <- score_ms2(scan_32, adenine, spec_entropy_params())
  expect_equal(score, e)
})