test_that("cluster_data works as expected", {
  mass_data_set <- readRDS(test_path("exttestdata", "test_mass_data_set.RDS"))
  distances <- dist_ms2(mass_data_set, 0.3, 2, gnps_params(0.5))
  results <- cluster_data(distances, mass_data_set, "opticlust")
  expect_true(nrow(results$abundance) == 6804)
  expect_true(ncol(results$abundance) == 3)
  expect_true(length(results) == 5)
})
