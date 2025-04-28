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

test_that("Test dist_shared works with bray", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)

  result <- dist_shared(communiy_object, 400, 10, "bray", 2)
  expect_true("data.frame" %in% class(result))
  expect_true(ncol(result) == 3)
  expect_true(nrow(result) == length(dat$samples)*length(dat$samples))
})

test_that("Test dist_shared works with jaccard", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)

  result <- dist_shared(communiy_object, 400, 10, "jaccard", 2)
  expect_true("data.frame" %in% class(result))
  expect_true(ncol(result) == 3)
  expect_true(nrow(result) == length(dat$samples)*length(dat$samples))
})

test_that("Test dist_shared works with hamming distance", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)

  result <- dist_shared(communiy_object, 400, 10, "hamming", 2)
  expect_true("data.frame" %in% class(result))
  expect_true(ncol(result) == 3)
  expect_true(nrow(result) == length(dat$samples)*length(dat$samples))
})

test_that("Test dist_shared works with soren index", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)

  result <- dist_shared(communiy_object, 400, 10, "soren", 2)
  expect_true("data.frame" %in% class(result))
  expect_true(ncol(result) == 3)
  expect_true(nrow(result) == length(dat$samples)*length(dat$samples))
})

test_that("Test dist_shared works with morisita horn index", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)

  result <- dist_shared(communiy_object, 400, 10, "morisita", 2)
  expect_true("data.frame" %in% class(result))
  expect_true(ncol(result) == 3)
  expect_true(nrow(result) == length(dat$samples)*length(dat$samples))
})

test_that("Test dist_shared works with thetayc(Yun and Clayton) distance", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)

  result <- dist_shared(communiy_object, 400, 10, "thetayc", 2)
  expect_true("data.frame" %in% class(result))
  expect_true(ncol(result) == 3)
  expect_true(nrow(result) == length(dat$samples)*length(dat$samples))
})


test_that("Test dist_shared errors when
          given the wrong community object", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  expect_error(dist_shared(result, 400, 10, "bray", 100))
})

test_that("Test dist_shared errors with wrong index", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  expect_error(dist_shared(result, 400, 10, "asad", 100))
})

test_that("Alpha summary returns the proper results for simpsons",{
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)
  alpha_sum <- alpha_summary(communiy_object, 400, 10, "simpson", 2)
  expect_true("data.frame" %in% class(alpha_sum))
  expect_true(ncol(alpha_sum) == length(dat$samples))
  expect_true(nrow(alpha_sum) == 1)
})

test_that("Alpha summary returns the proper results for shannon",{
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)
  alpha_sum <- alpha_summary(communiy_object, 400, 10, "shannon", 2)
  expect_true("data.frame" %in% class(alpha_sum))
  expect_true(ncol(alpha_sum) == length(dat$samples))
  expect_true(nrow(alpha_sum) == 1)
})

test_that("Alpha summary fails when given wrong input",{
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  results <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(results)
  expect_error(alpha_summary(results, 400, 10, "shannon", 2))
  expect_error(alpha_summary(communiy_object, 400, 10, "bray", 2))
})
