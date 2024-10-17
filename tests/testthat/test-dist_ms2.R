test_that("dist_ms2 works", {
  dir <- "exttestdata"
  file <- "demo_massdataset"
  dat <- readRDS(test_path(dir, file))

  first_10 <- dat@ms2_data[[1]]@variable_id[1:10]

  dat_sub <- dat %>%
    massdataset::activate_mass_dataset("variable_info") %>%
    massdataset::filter(variable_id %in% first_10)

  dist <- dist_ms2(dat_sub, cutoff = .3, precursor_thresh = 2,
                   gnps_params(frag_tolerance = 0.5))
  
  expect_s4_class(dist, "dgTMatrix")

  result <- data.frame(i = dist@i + 1,
                       j = dist@j + 1,
                       x = dist@x)
  
  expect_equal(nrow(result), 10)
})

test_that("dist_ms2_cpp works", {
  dir <- "exttestdata"
  file <- "demo_massdataset"
  dat <- readRDS(test_path(dir, file))

  first_10 <- dat@ms2_data[[1]]@variable_id[1:10]

  dat_sub <- dat %>%
    massdataset::activate_mass_dataset("variable_info") %>%
    massdataset::filter(variable_id %in% first_10)

  dist <- dist_ms2_cpp(dat_sub, 0.3, 2, gnps_params(0.5))
  expect_s3_class(dist, "data.frame")
  expect_equal(nrow(dist), 10)
})

test_that("dist_ms2_cpp gives the same output as dist_ms2", {
  dir <- "exttestdata"
  file <- "demo_massdataset"
  dat <- readRDS(test_path(dir, file))

  first_10 <- dat@ms2_data[[1]]@variable_id[1:10]

  dat_sub <- dat %>%
    massdataset::activate_mass_dataset("variable_info") %>%
    massdataset::filter(variable_id %in% first_10)

  dist <- dist_ms2(dat_sub, cutoff = .3, precursor_thresh = 2,
    gnps_params(frag_tolerance = 0.5))

  dist_r <- data.frame(i = dist@i + 1,
                       j = dist@j + 1,
                       dist = dist@x)
  
  dist_cpp <- dist_ms2_cpp(dat_sub, 0.3, 2, gnps_params(0.5))

  expect_true(all(dist_r == dist_cpp))
})



# library(testthat)
# dir <- "exttestdata"
# file <- "demo_massdataset"

# dat <- readRDS(test_path(dir, file))
# test_list <- list("pmz" = dat@ms2_data[[1]]@ms2_mz[1:10],
#                   "id" = dat@ms2_data[[1]]@variable_id[1:10],
#                   "spectra" = lapply(dat@ms2_data[[1]]@ms2_spectra[1:10],
#                                      as.data.frame))
# p <- gnps_params(frag_tolerance = 0.5)
# distMS2(test_list, p, 2, 0.7)
# benchmarking
# dir <- "exttest-data"
# file <- "demo_massdataset"
# dat <- readRDS(test_path(dir, file))

# first_10 <- dat@ms2_data[[1]]@variable_id[1:10]
# first_100 <- dat@ms2_data[[1]]@variable_id[1:100]
# first_500 <- dat@ms2_data[[1]]@variable_id[1:500]

# dat_sub_10 <- dat %>%
#     massdataset::activate_mass_dataset("variable_info") %>%
#     massdataset::filter(variable_id %in% first_10)

# dat_sub_100 <- dat %>%
#     massdataset::activate_mass_dataset("variable_info") %>%
#     massdataset::filter(variable_id %in% first_100)

# dat_sub_500 <- dat %>%
#     massdataset::activate_mass_dataset("variable_info") %>%
#     massdataset::filter(variable_id %in% first_500)
    
# microbenchmark::microbenchmark(dist_ms2(dat_sub_10, cutoff = .3, score_method = "gnps", frag_tolerance = 0.5),
#                                dist_ms2(dat_sub_100, cutoff = .3, score_method = "gnps", frag_tolerance = 0.5),
#                                dist_ms2(dat_sub_500, cutoff = .3, score_method = "gnps", frag_tolerance = 0.5),
#                                times = 2)


# profvis::profvis({dist_ms2(dat_sub_500, cutoff = .3, score_method = "gnps", precursor_thresh = 2, frag_tolerance = 0.5)})



