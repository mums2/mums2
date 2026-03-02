test_that("all the filters work as expected", {
  data <- import_all_data(peak_table = test_path("exttestdata",
                                                 "peak_table.csv"),
                          meta_data = test_path("exttestdata",
                                                "meta_data.csv"),
                          format = "Progenesis")
  current_peak_table <- get_peak_table(data)
  data_filtered <- data |>
    filter_peak_table(filter_mispicked_ions_params()) |>
    filter_peak_table(filter_cv_params(cv_threshold = 0.2)) |>
    filter_peak_table(filter_group_params(group_threshold = 0.1,
                                          "Blanks")) |>
    filter_peak_table(filter_insource_ions_params())
  filtered_peak_table <- get_peak_table(data_filtered)
  expect_false(nrow(current_peak_table) == nrow(filtered_peak_table))
  expect_true(nrow(current_peak_table) == 1303)
  expect_true(nrow(filtered_peak_table) == 590)
})

test_that("filter_peak_table will fail if given incorrect parameters", {
  data <- import_all_data(peak_table = test_path("exttestdata",
                                                 "peak_table.csv"),
                          meta_data = test_path("exttestdata",
                                                "meta_data.csv"),
                          format = "Progenesis")
  expect_error(filter_peak_table("a", filter_mispicked_ions_params()),
               "The mpactr object")
  expect_error(filter_peak_table(data, c()),
               "The params object")
})

test_that("filter_mispicked_ions parameters return correct data", {
  params <- filter_mispicked_ions_params()
  ls <- list(ringwin = 0, isowin = 0,
    trwin = 0, max_iso_shift = 0,
    merge_peaks = 0, merge_method = 0,
    copy_object = 0
  )
  expect_true(all(names(ls) == names(params)))
  expect_true(length(params) == 7)
})

test_that("filter_groups_parameters return correct data", {
  params <- filter_group_params(group_to_remove = "Blank")
  ls <- list(group_threshold = 0,
    group_to_remove = 0, remove_ions = 0,
    copy_object = 0
  )
  expect_true(all(names(ls) == names(params)))
  expect_true(length(params) == 4)
})

test_that("filter_cv_params return correct data", {
  params <- filter_cv_params(cv_threshold = 0.2)
  ls <- list(cv_threshold = 0,
             copy_object = FALSE)
  expect_true(all(names(ls) == names(params)))
  expect_true(length(params) == 2)
})

test_that("filter_insource_ions_params return correct data", {
  params <- filter_insource_ions_params(0.2)
  ls <- list(cluster_threshold = 0,
             copy_object = 0)
  expect_true(all(names(ls) == names(params)))
  expect_true(length(params) == 2)
})

