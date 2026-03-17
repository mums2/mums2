test_that("We can add ms2 data to our massdataset with mgf files", {
  data <- import_all_data(peak_table = test_path("exttestdata",
                                                 "peak_table.csv"),
                          meta_data = test_path("exttestdata",
                                                "meta_data.csv"),
                          format = "Progenesis")

  data_filtered <- data |>
    filter_peak_table(filter_mispicked_ions_params()) |>
    filter_peak_table(filter_cv_params(cv_threshold = 0.2)) |>
    filter_peak_table(filter_group_params(group_threshold = 0.1,
                                          "Blanks")) |>
    filter_peak_table(filter_insource_ions_params())

  mgf_file <-  test_path("exttestdata",
                         "12152023_Coculture_with_new_JC1.gnps.mgf")
  ms2_matches <- ms2_ms1_compare(mgf_file, data, 100, 150)
  expect_true(nrow(ms2_matches$ms2_matches) > 0)
})

test_that("We can add ms2 data to our massdataset with mzml files", {
  data <- import_all_data(peak_table = test_path("exttestdata",
                                                 "peak_table.csv"),
                          meta_data = test_path("exttestdata",
                                                "meta_data.csv"),
                          format = "Progenesis")

  data_filtered <- data |>
    filter_peak_table(filter_mispicked_ions_params()) |>
    filter_peak_table(filter_cv_params(cv_threshold = 0.2)) |>
    filter_peak_table(filter_group_params(group_threshold = 0.1,
                                          "Blanks")) |>
    filter_peak_table(filter_insource_ions_params())

  mzxml_files <- test_path("exttestdata", "threonine_i2_e35_pH_tree.mzXML")
  ms2_matches <- ms2_ms1_compare(mzxml_files, data, 100000, 150)
  expect_true(nrow(ms2_matches$ms2_matches) > 0)
})


test_that("ms2_ms2_compare fails if given incorrect parameters", {
  data <- import_all_data(peak_table = test_path("exttestdata",
                                                 "peak_table.csv"),
                          meta_data = test_path("exttestdata",
                                                "meta_data.csv"),
                          format = "Progenesis")

  data_filtered <- data |>
    filter_peak_table(filter_mispicked_ions_params()) |>
    filter_peak_table(filter_cv_params(cv_threshold = 0.2)) |>
    filter_peak_table(filter_group_params(group_threshold = 0.1,
                                          "Blanks")) |>
    filter_peak_table(filter_insource_ions_params())

  mzxml_files <- test_path("exttestdata", "threonine_i2_e35_pH_tree.mzXML")
  expect_error(ms2_ms1_compare(mzxml_files, "data", 100000, 150),
               "The mpactr object must be created using the")
  expect_error(ms2_ms1_compare(1234, data, 100000, 150),
               "ms2_files must be a character")
  expect_error(ms2_ms1_compare(mzxml_files, data, "100000", 150),
               "mz_tolerance must be numeric")
  expect_error(ms2_ms1_compare(mzxml_files, data, 100000, "150"),
               "rt_tolerance must be numeric")
})
