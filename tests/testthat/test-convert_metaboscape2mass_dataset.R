####################################################
####       test checks for metaboscape_ft       ####
####################################################

test_that("convert function aborts when a required sample_info
           column is missing", {
            meta_f <- "sample_metaboscape_metadata.csv"
            meta <- read.csv(test_path("exttestdata", meta_f))
            meta <- meta[, -3]

            file <- "sample_metaboscape.csv"
            df <- read.csv(test_path("exttestdata", file))

            expect_error(convert_metaboscape2mass_dataset(df, meta),
                         "are required columns that")
          })


test_that("convert function alerts about PEPMASS conversion", {
  meta_f <- "sample_metaboscape_metadata.csv"
  meta <- read.csv(test_path("exttestdata", meta_f))

  file <- "sample_metaboscape.csv"
  df <- read.csv(test_path("exttestdata", file))

  expect_message(convert_metaboscape2mass_dataset(df, meta),
                 "Compound mass was found as PEPMASS")
})

test_that("convert function aborts if ADDUCT column is not found", {
  meta_f <- "sample_metaboscape_metadata.csv"
  meta <- read.csv(test_path("exttestdata", meta_f))

  file <- "sample_metaboscape.csv"
  df <- read.csv(test_path("exttestdata", file))
  df <- df[, -5]

  expect_error(convert_metaboscape2mass_dataset(df, meta),
               "Compund mass was found as PEPMASS, but column ADDUCT")
})


####################################################
####  different class types for metaboscape_ft  ####
####################################################
test_that("a file path can be supplied", {
  meta_f <- "sample_metaboscape_metadata.csv"
  meta <- read.csv(test_path("exttestdata", meta_f))

  file <- "sample_metaboscape.csv"
  path <- test_path("exttestdata", file)

  out <- convert_metaboscape2mass_dataset(path, meta)

  expect_s4_class(out, "mass_dataset")
  expect_equal(ncol(out@expression_data), 6)
  expect_equal(nrow(out@sample_info), 6)
  expect_equal(nrow(out@variable_info), 10)
})

test_that("a data.frame can be supplied", {
  meta_f <- "sample_metaboscape_metadata.csv"
  meta <- read.csv(test_path("exttestdata", meta_f))

  file <- "sample_metaboscape.csv"
  df <- read.csv(test_path("exttestdata", file))

  out <- convert_metaboscape2mass_dataset(df, meta)

  expect_s4_class(out, "mass_dataset")
  expect_equal(ncol(out@expression_data), 6)
  expect_equal(nrow(out@sample_info), 6)
  expect_equal(nrow(out@variable_info), 10)
})


##################################################
####       test helper functions              ####
##################################################

test_that("get_variable_info returns tidymass formatted table", {
  meta_f <- "sample_metaboscape_metadata.csv"
  meta <- read.csv(test_path("exttestdata", meta_f))

  file <- "sample_metaboscape.csv"
  df <- read.csv(test_path("exttestdata", file))

  vi <- get_variable_info(df, meta)

  expect_true(all(c("variable_id", "rt", "mz") %in% colnames(vi)))
  expect_s3_class(vi, "data.frame")

})

test_that("mpactr_object properly converts to a massdataset", {
  peak_table_path <- mpactr::example_path("cultures_metaboscape_peaktable.csv")
  meta_data_path <- mpactr::example_path("cultures_metaboscape_metadata.csv" )
  mpact_data <- import_data(peak_table_path, meta_data_path, "Metaboscape")
  data_set <- convert_mpactr_object_to_mass_data_set(mpact_data)
  expect_s4_class(data_set, "mass_dataset")
})

test_that("conversion to massdataset to fail from wrong input", {
  peak_table_path <- mpactr::example_path("cultures_metaboscape_peaktable.csv")
  meta_data_path <- mpactr::example_path("cultures_metaboscape_metadata.csv" )
  mpact_data <- import_data(peak_table_path, meta_data_path, "Metaboscape")
  expect_error(convert_mpactr_object_to_mass_data_set(c()))
})
