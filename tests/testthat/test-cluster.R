test_that("cluster_data works as expected", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5), min_peaks = 0)
  results <- cluster_data(distances, dat,  0.3, "opticlust")
  expect_true(nrow(results$abundance) == 7119)
  expect_true(ncol(results$abundance) == 3)
  expect_true(length(results) == 5)
})

test_that("cluster_data error catching works as expected", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5), min_peaks = 0)
  expect_error(cluster_data("distances", dat,  0.3, "opticlust"),
               "distance_df must be an object created")
  expect_error(cluster_data(distances, "dat",  0.3, "opticlust"),
              "The mass_data object must be created")
  expect_error(cluster_data(distances, dat,  "0.3", "opticlust"),
              "cutoff should be a numeric value")
  
  zero_dist <- dist_ms2(dat, -1, 2, modified_cosine_params(0.5), min_peaks = 0)
  expect_error(cluster_data(zero_dist, dat,  0.3, "opticlust"),
              "distance_df must have more than 0 rows")
 
})

