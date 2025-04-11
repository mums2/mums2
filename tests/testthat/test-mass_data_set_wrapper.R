test_that("We can add ms2 data to our massdataset with mgf files", {
  data <- import_all_data(peak_table = test_path("exttestdata", "peak_table.csv"), 
                          meta_data = test_path("exttestdata", "meta_data.csv"), 
                          format = "Progenesis")

  data_filtered <- data |>
    filter_peak_table(filter_mispicked_ions_parameters()) |>
    filter_peak_table(filter_cv_parameters(cv_threshold = 0.2, cv_param = "mean")) |>
    filter_peak_table(filter_group_parameters(group_threshold = 0.1, "Blanks")) |>
    filter_peak_table(filter_insource_ions_parameters())

  mgf_files <- list.files(test_path("exttestdata"), pattern = ".mgf", full.names = TRUE)
  ms2_matches <- ms2_ms1_compare(mgf_files, data, 100, 150)
  expect_true(nrow(ms2_matches$ms2_matches) > 0)
})

test_that("We can add ms2 data to our massdataset with mzml files", {
  data <- import_all_data(peak_table = test_path("exttestdata", "peak_table.csv"), 
                          meta_data = test_path("exttestdata", "meta_data.csv"), 
                          format = "Progenesis")

  data_filtered <- data |>
    filter_peak_table(filter_mispicked_ions_parameters()) |>
    filter_peak_table(filter_cv_parameters(cv_threshold = 0.2, cv_param = "mean")) |>
    filter_peak_table(filter_group_parameters(group_threshold = 0.1, "Blanks")) |>
    filter_peak_table(filter_insource_ions_parameters())

  mgf_files <- list.files(test_path("exttestdata"), pattern = ".mgf", full.names = TRUE)
  ms2_matches <- ms2_ms1_compare(mgf_files, data, 100, 150)
  expect_true(nrow(ms2_matches$ms2_matches) > 0)
})
