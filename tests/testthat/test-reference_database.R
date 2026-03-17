test_that("We can get one of the reference by index from the database", {
  psu_msmls_data <- read_msp(test_path("exttestdata/database_data",
                                       "PSU-MSMLS.msp"))
  data <- get_reference_data(psu_msmls_data, 2)
  expect_equal(names(data), c("info", "spec"))
  expect_true(length(data$info$keys) == 15)
  expect_true(length(data$spec$mz) == 43)
  expect_error(get_reference_data(psu_msmls_data, ""),
               "Index has to be a numeric")
  expect_error(get_reference_data("", 1),
               "Ensure reference is the object")
})

test_that("get_reference_data index will fail if the size is greater than
          the length of the database", {
            psu_msmls_data <- read_msp(test_path("exttestdata/database_data",
                                                 "PSU-MSMLS.msp"))
            expect_error(get_reference_data(psu_msmls_data, 577),
                         "Index must be less than the size of the database")
          })


test_that("You can add another database file to the reference data", {
  path <- test_path("exttestdata/database_data", "PSU-MSMLS.msp")
  psu_msmls_data <- read_msp(path)
  psu_msmls_data2 <- read_msp(path)
  psu_msmls_data_combined <-
    combined_reference_database(psu_msmls_data, psu_msmls_data2)
  expect_equal(length(psu_msmls_data_combined), 1152)
  expect_error(combined_reference_database(psu_msmls_data, path),
               "Ensure reference is the object")
  expect_error(combined_reference_database("", psu_msmls_data_combined),
               "Ensure reference is the object")
})


test_that("We can add another reference to the database", {
  psu_msmls_data <- read_msp(test_path("exttestdata/database_data",
                                       "PSU-MSMLS.msp"))
  expect_output(print(psu_msmls_data),
                "You have 576 references in this object.")
})

test_that("Expect print to output a custom message", {
  psu_msmls_data <- read_msp(test_path("exttestdata/database_data",
                                       "PSU-MSMLS.msp"))
  expect_output(print(psu_msmls_data),
                "You have 576 references in this object.")
})
