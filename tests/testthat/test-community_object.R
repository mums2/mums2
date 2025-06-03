test_that("test that we can create a communiy object", {
  ms2_match_data <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  dist <- dist_ms2(ms2_match_data, 0.3, 2, gnps_params(0.5))
  final_cluster <- cluster_data(dist, ms2_match_data, 0.03, "opticlust")
  communiy_object <- create_community_matrix_object(final_cluster)
  expect_true("community_object" %in% class(communiy_object))
})

test_that("test that we can create a communiy object without clustering", {
  ms2_match_data <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  communiy_object <- create_community_matrix_object(ms2_match_data)
  expect_true("community_object" %in% class(communiy_object))
})

test_that("Test that get community object returns a community matrix", {
  ms2_match_data <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  dist <- dist_ms2(ms2_match_data, 0.3, 2, gnps_params(0.5))
  final_cluster <- cluster_data(dist, ms2_match_data, 0.03, "opticlust")
  communiy_object <- create_community_matrix_object(final_cluster)
  community_matrix <- get_community_matrix(communiy_object)
  expect_true("matrix" %in% class(community_matrix))
})
test_that("Printing a community object prints out the matrix", {
    ms2_match_data <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  dist <- dist_ms2(ms2_match_data, 0.3, 2, gnps_params(0.5))
  final_cluster <- cluster_data(dist, ms2_match_data, 0.03, "opticlust")
  communiy_object <- create_community_matrix_object(final_cluster)
  expect_output(print(communiy_object))
})

test_that("Expect get_community_object and matrix to error when given the wrong object", {
  count_table <- test_path("exttestdata", "final.count_table")
  expect_error(create_community_matrix_object(count_table))
  expect_error(get_community_matrix(count_table))
})

