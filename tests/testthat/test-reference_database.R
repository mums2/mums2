test_that("We can get one of the reference by index from the database", {
  psu_msmls_data <- read_msp(test_path("exttestdata/database_data",
                                       "PSU-MSMLS.msp"))
  data <- get_reference_data(psu_msmls_data, 1)
  expect_equal(names(data), c("info", "spec"))
  expect_true(length(data$info$keys) == 15)
  expect_true(length(data$spec$mz) == 43)
  expect_error(get_reference_data(psu_msmls_data, ""), "index has to be a numeric")
  expect_error(get_reference_data("", 1), "Ensure reference is the object")
})

test_that("You can add another database file to the reference data", {
  path <- test_path("exttestdata/database_data", "PSU-MSMLS.msp")
  psu_msmls_data <- read_msp(path)
  expect_equal(length(psu_msmls_data), 576)
  psu_msmls_data <- add_references(psu_msmls_data, path, "msp")
  expect_equal(length(psu_msmls_data), 1152)
  expect_error(add_references(psu_msmls_data, path, ""), "method has to be")
  expect_error(add_references("", "", ""), "Ensure reference is the object")
})


test_that("We can add another reference to the database", {
  psu_msmls_data <- read_msp(test_path("exttestdata/database_data",
                                      "PSU-MSMLS.msp"))
  expect_output(print(psu_msmls_data), "You have 576 references in this object.")
})

test_that("Expect print to output a custom message", {
  psu_msmls_data <- read_msp(test_path("exttestdata/database_data",
                                      "PSU-MSMLS.msp"))
  expect_output(print(psu_msmls_data), "You have 576 references in this object.")
})

