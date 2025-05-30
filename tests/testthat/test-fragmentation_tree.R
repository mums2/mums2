test_that("test that fragmentation tree makes proper predictions", {
  dat <- readRDS(test_path("exttestdata", "small_matched_data.RDS"))
  dat <- compute_molecular_formulas(dat)
  correct_results <- readRDS(test_path("exttestdata", "prediction_molecular_formula_results.RDS"))
  expect_equal(dat$predicted_molecular_formulas, correct_results)
})
