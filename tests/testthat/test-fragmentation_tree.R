test_that("test that fragmentation tree makes proper predictions", {
  dat <- readRDS(test_path("exttestdata", "small_matched_data.RDS"))
  dat <- compute_molecular_formulas(dat)
  correct_results <- readRDS(test_path("exttestdata", "prediction_molecular_formula_results.RDS"))
  expect_equal(dat$predicted_molecular_formulas, correct_results)
})

test_that("Warning messages appear and return empty when there are no parent decompositions", {
  dat <- readRDS(test_path("exttestdata", "small_matched_data.RDS"))
  dat$ms2_matches$mz[[1]] <- 70
  dat$ms2_matches <- dat$ms2_matches[1,]
  dat$peak_data <- dat$peak_data[1]
  expect_warning(compute_molecular_formulas(dat), "No parent decompositions")
})

test_that("Will return the first candidate if there is only one parent decomposition", {
  dat <- readRDS(test_path("exttestdata", "small_matched_data.RDS"))
  dat$ms2_matches$mz[[1]] <- 99
  dat$ms2_matches <- dat$ms2_matches[1,]
  dat$peak_data <- dat$peak_data[1]
  expect_warning(compute_molecular_formulas(dat, 1), "No valid parent formulas")
})

test_that("Will return the first valid index if there is only one", {
  dat <- readRDS(test_path("exttestdata", "small_matched_data.RDS"))
  dat$ms2_matches$mz[[1]] <- 99
  dat$ms2_matches <- dat$ms2_matches[1,]
  dat$peak_data <- dat$peak_data[1]
  result <- compute_molecular_formulas(dat, 2)
  Rdisop::decomposeMass(99)$formula[[2]]
  expect_equal(result$predicted_molecular_formulas[[1]],  Rdisop::decomposeMass(99)$formula[[2]])
})


test_that("Returns the first candidate if there are no children decompositions", {
  dat <- readRDS(test_path("exttestdata", "small_matched_data.RDS"))
  dat$ms2_matches$mz[1] <- 150
  dat$ms2_matches <- dat$ms2_matches[1:2,]
  dat$peak_data <- dat$peak_data[1]
  dat$peak_data[[1]]$mz <- 70
  dat$peak_data[[1]]$intensity <- 10
  result <- compute_molecular_formulas(dat)
  expected_result <- Rdisop::decomposeMass(150)$formula[[1]]
  expect_equal(result$predicted_molecular_formulas[[1]], expected_result)
})


test_that("Returns replaces valid indexes with invalid indexes if there are no valid indexes", {
  dat <- readRDS(test_path("exttestdata", "small_matched_data.RDS"))
  dat$ms2_matches$mz[1] <- 130
  dat$ms2_matches <- dat$ms2_matches[1:2,]
  dat$peak_data <- dat$peak_data[1]
  dat$peak_data[[1]]$mz <- 70
  dat$peak_data[[1]]$intensity <- 10
  expect_warning(compute_molecular_formulas(dat), "No valid parent formulas")
})


