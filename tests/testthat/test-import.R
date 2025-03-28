test_that("Ensure we can import data properly", {
  data <- import_all_data(peak_table = test_path("exttestdata", "peak_table.csv"), 
  meta_data = test_path("exttestdata", "meta_data.csv"), 
  format = "Progenesis")
  
  expect_true("filter_pactr" %in% class(data))
  expect_error(import_all_data(peak_table = "", 
  meta_data = test_path("exttestdata", "meta_data.csv"), 
  format = "Progenesis"))

  expect_error(import_all_data(peak_table = test_path("exttestdata", "peak_table.csv"), 
  meta_data = data.table(), 
  format = "Progenesis"))
})
