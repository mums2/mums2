test_that("I can get all slots within the mass_data object", {
  dat <- readRDS(test_path("exttestdata", "small_matched_data.RDS"))

  expect_equal(dat$ms2_matches, get_ms2_matches(dat))
  expect_equal(dat$ms1_data, get_ms1_data(dat))
  expect_equal(dat$peak_data, get_ms2_peaks_data(dat))
  expect_equal(dat$samples, get_samples(dat))
  expect_equal(dat$predicted_molecular_formulas,
               get_molecular_formula_preds(dat))
})

test_that("Getter functions fail if not given the proper parameter", {
  dat <- readRDS(test_path("exttestdata", "small_matched_data.RDS"))
  dat$predicted_molecular_formulas <- NULL
  expect_error(get_molecular_formula_preds(dat),
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
  expect_error(get_molecular_formula_preds(""),
               "mass_data must be generated from the")
})
