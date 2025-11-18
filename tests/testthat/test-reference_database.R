test_that("We can get one of the reference by index from the database", {
  psu_msmls_data <- read_msp(test_path("exttestdata/database_data",
                                       "PSU-MSMLS.msp"))
  data <- get_reference_data(psu_msmls_data, 1)
  expect_equal(names(data), c("info", "spec"))
  expect_true(length(data$info$keys) == 15)
  expect_true(length(data$spec$mz) == 43)
})

test_that("Expect print to output a custom message", {
  psu_msmls_data <- read_msp(test_path("exttestdata/database_data",
                                      "PSU-MSMLS.msp"))
  expect_output(print(psu_msmls_data), "You have 576 references in this object.")
})