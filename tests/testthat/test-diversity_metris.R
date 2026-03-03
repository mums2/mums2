test_that("Test dist_shared works with bray", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  community_object <- create_community_matrix_object(result)

  result <- dist_shared(community_object, 400, 10, "bray", TRUE,
                        iterations = 2)
  expect_true("dist" %in% class(result))
  expect_true(length(result) == 210)
})

test_that("Test dist_shared works with without subsample = FALSE", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  community_object <- create_community_matrix_object(result)

  result <- dist_shared(community_object, 400, 10, "bray", FALSE,
                        iterations = 2)
  expect_true("dist" %in% class(result))
  expect_true(length(result) == 210)
})


test_that("Test dist_shared works with jaccard", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  community_object <- create_community_matrix_object(result)

  result <- dist_shared(community_object, 400, 10, "jaccard", TRUE,
                        iterations = 2)
  expect_true("dist" %in% class(result))
  expect_true(length(result) == 210)
})

test_that("Test dist_shared works with hamming distance", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  community_object <- create_community_matrix_object(result)

  result <- dist_shared(community_object, 400, 10, "hamming", TRUE,
                        iterations = 2)
  expect_true("dist" %in% class(result))
  expect_true(length(result) == 210)
})

test_that("Test dist_shared works with soren index", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  community_object <- create_community_matrix_object(result)

  result <- dist_shared(community_object, 400, 10, "soren", TRUE,
                        iterations = 2)
  expect_true("dist" %in% class(result))
  expect_true(length(result) == 210)
})

test_that("Test dist_shared works with morisita horn index", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  community_object <- create_community_matrix_object(result)

  result <- dist_shared(community_object, 400, 10, "morisita", TRUE,
                        iterations = 2)
  expect_true("dist" %in% class(result))
  expect_true(length(result) == 210)
})

test_that("Test dist_shared works with thetayc(Yun and Clayton) distance", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  community_object <- create_community_matrix_object(result)

  result <- dist_shared(community_object, 400, 10, "thetayc", TRUE,
                        iterations = 2)
  expect_true("dist" %in% class(result))
  expect_true(length(result) == 210)
})


test_that("Test dist_shared errors when
          given the wrong community object", {
            dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
            distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
            result <- cluster_data(distances, dat,  0.3, "opticlust")
            expect_error(dist_shared(result, 400, 10, "bray", TRUE,
                                     iterations = 100))
          })

test_that("Test dist_shared errors with wrong object", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  expect_error(dist_shared(result, 400, 10, "asad", TRUE,
                           iterations = 100))
})

test_that("Test dist_shared errors with wrong index", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  community_object <- create_community_matrix_object(result)
  expect_error(dist_shared(community_object, 400, 10, "asad", TRUE,
                           iterations = 100))
})

test_that("Test dist_shared errors when given wrong parameters", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  community_object <- create_community_matrix_object(result)
  expect_error(dist_shared(community_object, "400", 10, "bray", TRUE,
                           iterations = 100),
               "size")
  expect_error(dist_shared(community_object, 400, "10", "bray", TRUE,
                           iterations = 100),
              "threshold")
  expect_error(dist_shared(community_object, 400, 10, "bray", "TRUE",
                           iterations = 100),
              "subsample")
  expect_error(dist_shared(community_object, 400, 10, "bray", TRUE,
                           iterations = "100"),
              "iterations")
  expect_error(dist_shared(community_object, 400, 10, "bray", TRUE,
                           iterations = 100, number_of_threads = "1"),
              "number_of_threads")
  expect_error(dist_shared(community_object, 400, 10, "bray", TRUE,
                           iterations = 100, seed = "1"),
              "seed")

})

test_that("Alpha summary returns the proper results for simpsons", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  community_object <- create_community_matrix_object(result)
  alpha_sum <- alpha_summary(community_object, 400, 10, "simpson", TRUE,
                             iterations = 2)
  expect_true("matrix" %in% class(alpha_sum))
  expect_true(ncol(alpha_sum) == length(dat$samples))
  expect_true(nrow(alpha_sum) == 1)
})

test_that("Alpha summary returns the proper results for shannon", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  community_object <- create_community_matrix_object(result)
  alpha_sum <- alpha_summary(community_object, 400, 10, "shannon", TRUE,
                             iterations = 2)
  expect_true("matrix" %in% class(alpha_sum))
  expect_true(ncol(alpha_sum) == length(dat$samples))
  expect_true(nrow(alpha_sum) == 1)
})

test_that("Alpha summary works when subsample = FALSE", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
  result <- cluster_data(distances, dat,  0.3, "opticlust")
  community_object <- create_community_matrix_object(result)
  alpha_sum <- alpha_summary(community_object, 400, 10, "simpson", FALSE,
                             iterations = 2)
  expect_true("matrix" %in% class(alpha_sum))
  expect_true(ncol(alpha_sum) == length(dat$samples))
  expect_true(nrow(alpha_sum) == 1)
})


test_that("Alpha summary fails when given wrong input", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
  results <- cluster_data(distances, dat,  0.3, "opticlust")
  community_object <- create_community_matrix_object(results)
  expect_error(alpha_summary(results, 400, 10, "shannon", TRUE, 2))
  expect_error(alpha_summary(community_object, 400, 10, "bray", TRUE,
                             iterations = 2))
})

test_that("Test Alpha summary errors when given wrong parameters", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
  results <- cluster_data(distances, dat,  0.3, "opticlust")
  community_object <- create_community_matrix_object(results)
  expect_error(alpha_summary(community_object, "400", 10, "shannon", TRUE,
                           iterations = 100),
               "size")
  expect_error(alpha_summary(community_object, 400, "10", "shannon", TRUE,
                           iterations = 100),
              "threshold")
  expect_error(alpha_summary(community_object, 400, 10, "shannon", "TRUE",
                           iterations = 100),
              "subsample")
  expect_error(alpha_summary(community_object, 400, 10, "shannon", TRUE,
                           iterations = "100"),
              "iterations")
  expect_error(alpha_summary(community_object, 400, 10, "shannon", TRUE,
                           iterations = 100, number_of_threads = "1"),
              "number_of_threads")
  expect_error(alpha_summary(community_object, 400, 10, "shannon", TRUE,
                           iterations = 100, seed = "1"),
              "seed")

})
