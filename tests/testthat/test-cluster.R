test_that("cluster_data works as expected", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  results <- cluster_data(distances, dat,  0.3, "opticlust")
  expect_true(nrow(results$abundance) == 7119)
  expect_true(ncol(results$abundance) == 3)
  expect_true(length(results) == 5)
})
