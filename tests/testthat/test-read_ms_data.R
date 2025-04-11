test_that("read_mgf will read an mgf data properly", {
  mgf_files <- list.files(test_path("exttestdata"), pattern = ".mgf", full.names = TRUE)
  mgf_data <- read_mgf(mgf_files)
  expect_true(length(mgf_data) == 2)
  expect_true(length(mgf_data$peak_data[[1]]) == nrow(mgf_data$mass_spec_data)) 
})

test_that("read_msp will read an msp data properly", {
  psu_msmls_data <- read_msp(test_path("exttestdata/database_data", "PSU-MSMLS.msp"))
  expect_true(length(psu_msmls_data) == 576)
  expect_true(class(psu_msmls_data) %in% "list")
})

