test_that("read_mgf will read an mgf data properly", {
  mgf_files <- list.files(test_path("exttestdata"), pattern = ".mgf", full.names = TRUE)
  mgf_data <- read_mgf(mgf_files)
  expect_true(length(mgf_data) == 2)
  expect_true(length(mgf_data$peak_data[[1]]) == nrow(mgf_data$mass_spec_data)) 
})

test_that("read_mgf will fail if file has the wrong extension", {
  expect_error(read_mgf(test_path("exttestdata", "matched_data.RDS")))
})

test_that("read_mgf will fail if the file does not exist", {
  expect_error(read_mgf(".mgf"))
})

test_that("read_mzml_mzxml will read an mgf data properly", {
  path <- test_path("exttestdata", "threonine_i2_e35_pH_tree.mzXML")
  mzml_data <- read_mzml_mzxml(path)
  expect_true(length(mzml_data) == 2)
  expect_true(length(mzml_data$peak_data[[1]]) == nrow(mzml_data$mass_spec_data)) 
})

test_that("read_mzml_mzxml will fail if file has the wrong extension", {
  expect_error(read_mzml_mzxml(test_path("exttestdata", "matched_data.RDS")))
})

test_that("read_mzml_mzxml will fail if the file does not exist", {
  expect_error(read_mzml_mzxml(".mzml"))
  expect_error(read_mzml_mzxml(".mzxml"))
})

test_that("read_msp will read an msp data properly", {
  psu_msmls_data <- read_msp(test_path("exttestdata/database_data", "PSU-MSMLS.msp"))
  expect_true(length(psu_msmls_data) == 576)
  expect_true(class(psu_msmls_data) %in% "list")
})

test_that("read_msp will fail if file has the wrong extension", {
  expect_error(read_msp(test_path("exttestdata", "matched_data.RDS")))
})

