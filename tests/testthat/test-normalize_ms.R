test_that("relative_abundance works", {
  raw_values <- c(0.00, 1.00, 0.00, 5.00, 2.50, 1.50)

  norm_values <- relative_abundance(raw_values)
  expect_equal(norm_values, c(0, 0.1, 0, 0.5, 0.25, 0.15))
})

test_that("normalize_feature_table produces error if the expression_data is not activated", {
  dir <- "exttestdata"
  file <- "demo_massdataset"
  dat <- readRDS(test_path(dir, file))

  expect_error(normalize_ms(dat, "rel_abund"),
               "activate your object using activate_mass_dataset")
})

test_that("normalize_feature_table produces error if NA are in expression_data", {
  dir <- "exttestdata"
  file <- "demo_massdataset"
  dat <- readRDS(test_path(dir, file))

  expect_error(dat %>%
                massdataset::activate_mass_dataset(what = "expression_data") %>% 
                normalize_ms(method = "rel_abund"),
              "Expression data contains NA values")
})

test_that("normalize_feature_table works on mass_dataset with method = rel_abund", {

  dir <- "exttestdata"
  file <- "demo_massdataset"
  dat <- readRDS(test_path(dir, file))

  dat@expression_data[is.na(dat@expression_data)] <- 0
  dat_norm <- dat %>%
    massdataset::activate_mass_dataset("expression_data") %>%
    normalize_ms("rel_abund")

  sum_103 <- sum(dat@expression_data[, "sample_103"])
  expected <- dat@expression_data[, "sample_103"] / sum_103
  
  expect_s4_class(dat_norm, "mass_dataset")
  expect_equal(round(dat_norm@expression_data[1:14, "sample_103"], 8), round(expected[1:14], 8))
  expect_equal(dat_norm@process_info[[length(dat_norm@process_info)]]@function_name, "normalize_ms()")

})

