test_that("Ensure we can import data properly", {
  data <- import_all_data(peak_table = test_path("exttestdata",
                                                 "peak_table.csv"),
                          meta_data = test_path("exttestdata",
                                                "meta_data.csv"),
                          format = "Progenesis")

  expect_true("filter_pactr" %in% class(data))
  expect_error(import_all_data(peak_table = "",
                               meta_data = test_path("exttestdata",
                                                     "meta_data.csv"),
                               format = "Progenesis"))

  expect_error(import_all_data(peak_table = test_path("exttestdata",
                                                      "peak_table.csv"),
                               meta_data = data.table(),
                               format = "Progenesis"))
})

test_that("Ensure data is properly converted to utf-8", {
  data <- mpactr::import_data(peak_table = test_path("exttestdata",
                                                     "peak_table.csv"),
                              meta_data = test_path("exttestdata",
                                                    "meta_data.csv"),
                              format = "Progenesis")
  dat <- get_peak_table(data)
  corrupted_feature_name <- "779_[Dâ€?Asp3,Dha7]MCâ€?FR_1000.50475 Da 187.19 s"
  dat$Compound[[1]] <-  corrupted_feature_name
  data$mpactr_data$set_peak_table(dat)
  data <- format_data_to_uft8_and_remove_commas(data)
  fixed_dat <- get_peak_table(data)
  fixed_feature <- fixed_dat$Compound[[1]]
  expect_true(fixed_feature != corrupted_feature_name)
})

test_that("Ensure we can change rt to RTINMINUTES or RTINSECONDS", {
  data <- import_all_data(peak_table = test_path("exttestdata",
                                                 "peak_table.csv"),
                          meta_data = test_path("exttestdata",
                                                "meta_data.csv"),
                          format = "Progenesis")
  data <- change_rt_to_seconds_or_minutes(data, "seconds")
  expect_true("RTINSECONDS" %in% colnames(get_peak_table(data)))
  data <- change_rt_to_seconds_or_minutes(data, "minutes")
  expect_true("RTINMINUTES" %in% colnames(get_peak_table(data)))
})

test_that("Ensure we change_rt_to_seconds_or_minutes fails if
          given the wrong object", {
            expect_error(change_rt_to_seconds_or_minutes("a", "seconds"))
          })

test_that("Ensure we change_rt_to_seconds_or_minutes fails if the data frame
          doesn't have correct columns", {
            data <- import_all_data(peak_table = test_path("exttestdata",
                                                           "peak_table.csv"),
                                    meta_data = test_path("exttestdata",
                                                          "meta_data.csv"),
                                    format = "Progenesis")
            data$mpactr_data$set_peak_table(data.frame())
            expect_error(change_rt_to_seconds_or_minutes(data, "seconds"))
            expect_error(change_rt_to_seconds_or_minutes(data, "minutes"))
          })
