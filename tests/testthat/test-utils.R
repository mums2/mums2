### get_peaks_data

test_that("get_peaks_data works", {
  dir <- "exttestdata"
  file <- "demo_massdataset"
  dat <- readRDS(test_path(dir, file))

  peaks <- get_peaks_data(dat, 1)

  e_mz <- dat@ms2_data[[1]]@ms2_mz[1]
  e_peaks <- dat@ms2_data[[1]]@ms2_spectra[[1]]

  expect_equal(length(peaks), 2)
  expect_equal(peaks$precursor_mz, e_mz)
  expect_true(all(peaks$spectra == e_peaks))
})

### create_peaks_data

test_that("create_peaks_data works", {
  file <- example("PSUMSMLS_Adenine.csv")
  peaks <- create_peaks_data(read.csv(file),
                             136.0620)

  expect_equal(names(peaks), c("precursor_mz", "spectra"))
  expect_s3_class(peaks, "peaks_data")
})
