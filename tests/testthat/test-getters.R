test_that("I can get all slots within the mass_data object", {
  dat <- readRDS(test_path("exttestdata", "small_matched_data.RDS"))

  expect_equal(dat$ms2_matches, get_ms2_matches(dat))
  expect_equal(dat$ms1_data, get_ms1_data(dat))
  expect_equal(dat$peak_data, get_ms2_peaks_data(dat))
  expect_equal(dat$samples, get_samples(dat))
  expect_equal(dat$predicted_molecular_formulas,
               get_molecular_formula_predictions(dat))
  
  data <- import_all_data(peak_table = test_path("exttestdata",
                                                 "peak_table.csv"),
                          metadata = test_path("exttestdata",
                                                "metadata.csv"),
                          format = "Progenesis")
  expect_equal(get_peak_table(data), get_peaks_table_mpactr(data))
  expect_equal(get_metadata(data), get_metadata_mpactr(data))
})

test_that("Getter functions fail if not given the proper parameter", {
  dat <- readRDS(test_path("exttestdata", "small_matched_data.RDS"))
  dat$predicted_molecular_formulas <- NULL
  expect_error(get_molecular_formula_predictions(dat),
               "No data found, make sure you run the")
  dat <- readRDS(test_path("exttestdata", "small_matched_data.RDS"))

  expect_error(get_ms2_matches(""),
               "mass_data must be generated from the")
  expect_error(get_ms1_data(""),
               "mass_data must be generated from the")
  expect_error(get_ms2_peaks_data(""),
               "mass_data must be generated from the")
  expect_error(get_samples(""),
               "mass_data must be generated from the")
  expect_error(get_molecular_formula_predictions(""),
               "mass_data must be generated from the")

  expect_error(get_peaks_table_mpactr(""),
               "mpactr object must be generated from the")
  expect_error(get_metadata_mpactr(""),
               "mpactr object must be generated from the")
})

