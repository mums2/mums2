test_that("Return error if diverstiy metric is not valid", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)
  expect_error(diversity(communiy_object, "no"))
  expect_error(diversity(c(), "shannon"))
})

test_that("Diversity returns a diversity value for every sample", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)
  samples <- rownames(communiy_object)
  diversity_result <- diversity(communiy_object, "shannon")
  expect_true(all(colnames(diversity_result) == samples))
})

test_that("Diversity metric shannon works", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)
  samples <- rownames(communiy_object)
  diversity_result <- diversity(communiy_object, "shannon")
  expect_true("matrix" %in% class(diversity_result))
})


test_that("Diversity metric simpson works", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)
  samples <- rownames(communiy_object)
  diversity_result <- diversity(communiy_object, "simpson")
  expect_true("matrix" %in% class(diversity_result))
})


test_that("Diversity metric bray works", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)
  samples <- rownames(communiy_object)
  diversity_result <- diversity(communiy_object, "bray")
  expect_true("matrix" %in% class(diversity_result))
  expect_true(ncol(diversity_result) == nrow(diversity_result))
})


test_that("Diversity errors when giving the wrong index", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)
  expect_error(diversity(communiy_object, "a"))
})

test_that("Test average_subsampled_dissimilarity works", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)

  result <- averaged_subsampled_dissimilarity(communiy_object, 400, 10, "bray", 100)
  expect_true("data.frame" %in% class(result))
  expect_true(ncol(result) == 3)
  expect_true(nrow(result) == length(dat$samples)*length(dat$samples))
})

test_that("Test average_subsampled_dissimilarity errors when
          given the wrong community object", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  expect_error(averaged_subsampled_dissimilarity(result, 400, 10, "bray", 100))
})
