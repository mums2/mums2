test_that("Return error if diverstiy metric is not valid", {
  count_table <- test_path("exttestdata", "final.count_table")
  distances <- test_path("exttestdata", "final.dist")
  final_count <- read_count(count_table)
  final_dist <- read_dist(distances, final_count, 0.03, F)
  final_cluster <- cluster(final_dist, 0.03, "opticlust")
  m <- create_community_matrix_object(final_cluster)
  expect_error(diversity(m, "no"))
  expect_error(diversity(c(), "shannon"))
})

test_that("Diversity returns a diversity value for every sample", {
  count_table <- test_path("exttestdata", "final.count_table")
  distances <- test_path("exttestdata", "final.dist")
  final_count <- read_count(count_table)
  final_dist <- read_dist(distances, final_count, 0.03, F)
  final_cluster <- cluster(final_dist, 0.03, "opticlust")
  m <- create_community_matrix_object(final_cluster)
  samples <- rownames(m)
  div <- diversity(m, "shannon")
  expect_true(all(colnames(div) == samples))
})

test_that("Diversity metric shannon works", {
  count_table <- test_path("exttestdata", "final.count_table")
  distances <- test_path("exttestdata", "final.dist")
  final_count <- read_count(count_table)
  final_dist <- read_dist(distances, final_count, 0.03, F)
  final_cluster <- cluster(final_dist, 0.03, "opticlust")
  m <- create_community_matrix_object(final_cluster)
  samples <- rownames(m)
  div <- diversity(m, "shannon")
  expect_true("matrix" %in% class(div))
  # need to save the results, will wait for actually data
})


test_that("Diversity metric simpson works", {
  count_table <- test_path("exttestdata", "final.count_table")
  distances <- test_path("exttestdata", "final.dist")
  final_count <- read_count(count_table)
  final_dist <- read_dist(distances, final_count, 0.03, F)
  final_cluster <- cluster(final_dist, 0.03, "opticlust")
  m <- create_community_matrix_object(final_cluster)
  samples <- rownames(m)
  div <- diversity(m, "simpson")
  expect_true("matrix" %in% class(div))
  # need to save the results, will wait for actually data
})


test_that("Diversity metric bray works", {
  count_table <- test_path("exttestdata", "final.count_table")
  distances <- test_path("exttestdata", "final.dist")
  final_count <- read_count(count_table)
  final_dist <- read_dist(distances, final_count, 0.03, F)
  final_cluster <- cluster(final_dist, 0.03, "opticlust")
  m <- create_community_matrix_object(final_cluster)
  samples <- rownames(m)
  div <- diversity(m, "bray")
  expect_true("matrix" %in% class(div))
  expect_true(ncol(div) == nrow(div))
  # need to save the results, will wait for actually data
})
