test_that("read_mgf will read an mgf data properly", {

  mgf_files <- test_path("exttestdata",
                         "12152023_Coculture_with_new_JC1.gnps.mgf")
  mgf_data <- read_mgf(mgf_files)
  expect_true(length(mgf_data) == 2)
  expect_true(length(mgf_data$peak_data[[1]]) == nrow(mgf_data$mass_spec_data))
})

test_that("read_mgf will fail if file has the wrong extension", {
  expect_error(read_mgf(test_path("exttestdata", "matched_data.RDS")))
})

test_that("read_mgf will fail if the file does not exist", {
  expect_error(read_mgf(".mgf"))
})

test_that("read_mzml_mzxml will read an mgf data properly", {
  path <- test_path("exttestdata", "threonine_i2_e35_pH_tree.mzXML")
  mzml_data <- read_mzml_mzxml(path)
  expect_true(length(mzml_data) == 2)
  expect_true(length(mzml_data$peak_data[[1]])
              == nrow(mzml_data$mass_spec_data))
})

test_that("read_mzml_mzxml will fail if file has the wrong extension", {
  expect_error(read_mzml_mzxml(test_path("exttestdata", "matched_data.RDS")))
})

test_that("read_mzml_mzxml will fail if the file does not exist", {
  expect_error(read_mzml_mzxml(".mzml"))
  expect_error(read_mzml_mzxml(".mzxml"))
})

test_that("read_msp will read an msp data properly", {
  psu_msmls_data <- read_msp(test_path("exttestdata/database_data",
                                       "PSU-MSMLS.msp"))
  expect_true(length(psu_msmls_data) == 576)
  expect_true(class(psu_msmls_data) %in% "reference_database")
})

test_that("read_msp will fail if file has the wrong extension", {
  expect_error(read_msp(test_path("exttestdata", "matched_data.RDS")))
})

test_that("read_hmdb will properly read hmdb data", {
  hmdb_db <- read_hmdb(test_path("exttestdata", "sweat_metabo_small.xml"),
                 test_path("exttestdata", "Spectra"))
  expect_true(length(hmdb_db) == 2)
})

test_that("read_hmdb will fail if given incorrect parameters", {
  expect_error(read_hmdb("", test_path("exttestdata", "Spectra")),
  "hmdb file does not exist")

  expect_error(read_hmdb(test_path("exttestdata", "sweat_metabo_small.xml"),
               "ms2_folder does not exist"))
})

test_that("read_hmdb will default prescursor mz to NA if not given", {
  hmdb_db <- read_hmdb(test_path("exttestdata", "sweat_metabo_small_no_mass.xml"),
                test_path("exttestdata", "Spectra"))

  node <- get_reference_data(hmdb_db, 1)
  expect_true(node$info$values[which(node$info$keys == "precursormz")] == "NA") 
})

