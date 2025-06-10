test_that("Test dist_shared works with bray", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)

  result <- dist_shared(communiy_object, 400, 10, "bray", T, 2)
  expect_true("dist" %in% class(result))
  expect_true(length(result) == 210)
})

test_that("Test dist_shared works with without subsample = F", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)

  result <- dist_shared(communiy_object, 400, 10, "bray", F, 2)
  expect_true("dist" %in% class(result))
  expect_true(length(result) == 210)
})


test_that("Test dist_shared works with jaccard", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)

  result <- dist_shared(communiy_object, 400, 10, "jaccard",T , 2)
  expect_true("dist" %in% class(result))
  expect_true(length(result) == 210)
})

test_that("Test dist_shared works with hamming distance", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)

  result <- dist_shared(communiy_object, 400, 10, "hamming", T, 2)
  expect_true("dist" %in% class(result))
  expect_true(length(result) == 210)
})

test_that("Test dist_shared works with soren index", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)

  result <- dist_shared(communiy_object, 400, 10, "soren", T, 2)
  expect_true("dist" %in% class(result))
  expect_true(length(result) == 210)
})

test_that("Test dist_shared works with morisita horn index", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)

  result <- dist_shared(communiy_object, 400, 10, "morisita", T, 2)
  expect_true("dist" %in% class(result))
  expect_true(length(result) == 210)
})

test_that("Test dist_shared works with thetayc(Yun and Clayton) distance", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)

  result <- dist_shared(communiy_object, 400, 10, "thetayc", T, 2)
  expect_true("dist" %in% class(result))
  expect_true(length(result) == 210)
})


test_that("Test dist_shared errors when
          given the wrong community object", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  expect_error(dist_shared(result, 400, 10, "bray", T, 100))
})

test_that("Test dist_shared errors with wrong object", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  expect_error(dist_shared(result, 400, 10, "asad", T, 100))
})

test_that("Test dist_shared errors with wrong index", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)
  expect_error(dist_shared(communiy_object, 400, 10, "asad", T, 100))
})

test_that("Alpha summary returns the proper results for simpsons",{
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)
  alpha_sum <- alpha_summary(communiy_object, 400, 10, "simpson", T, 2)
expect_true("matrix" %in% class(alpha_sum))
  expect_true(ncol(alpha_sum) == length(dat$samples))
  expect_true(nrow(alpha_sum) == 1)
})

test_that("Alpha summary returns the proper results for shannon",{
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)
  alpha_sum <- alpha_summary(communiy_object, 400, 10, "shannon", T, 2)
expect_true("matrix" %in% class(alpha_sum))
  expect_true(ncol(alpha_sum) == length(dat$samples))
  expect_true(nrow(alpha_sum) == 1)
})

test_that("Alpha summary works when subsample = F",{
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(result)
  alpha_sum <- alpha_summary(communiy_object, 400, 10, "simpson", F, 2)
  expect_true("matrix" %in% class(alpha_sum))
  expect_true(ncol(alpha_sum) == length(dat$samples))
  expect_true(nrow(alpha_sum) == 1)
})


test_that("Alpha summary fails when given wrong input",{
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  results <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(results)
  expect_error(alpha_summary(results, 400, 10, "shannon", T, 2))
  expect_error(alpha_summary(communiy_object, 400, 10, "bray", T, 2))
})
