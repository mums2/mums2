test_that("We can add ms2 data to our massdataset", {
  data <- import_all_data(peak_table = test_path("exttestdata", "peak_table.csv"), 
                          meta_data = test_path("exttestdata", "meta_data.csv"), 
                          format = "Progenesis")

  data_filtered <- data |>
    filter_peak_table(filter_mispicked_ions_parameters()) |>
    filter_peak_table(filter_cv_parameters(cv_threshold = 0.2, cv_param = "mean")) |>
    filter_peak_table(filter_group_parameters(group_threshold = 0.1, "Blanks")) |>
    filter_peak_table(filter_insource_ions_parameters())

  mass_data_set <- convert_mpactr_object_to_mass_data_set(data_filtered)
  
  mass_data_set <- add_ms2(mass_data_set = mass_data_set, 
    path = test_path("exttestdata"), 
    polarity = "negative", column = "hilic",
    ms1.ms2.match.mz.tol = 100,
    ms1.ms2.match.rt.tol = 150
  )
  expect_true(length(mass_data_set@ms2_data) > 0)
})
