test_that("dist_ms2_cpp works", {
  dir <- "exttestdata"
  file <- "demo_massdataset"
  dat <- readRDS(test_path(dir, file))

  first_10 <- dat@ms2_data[[1]]@variable_id[1:10]

  dat_sub <- dat %>%
    massdataset::activate_mass_dataset("variable_info") %>%
    massdataset::filter(variable_id %in% first_10)

  dist <- dist_ms2(dat_sub, 0.3, 2, gnps_params(0.5))
  expect_s3_class(dist, "data.frame")
  expect_equal(nrow(dist), 10)
})




