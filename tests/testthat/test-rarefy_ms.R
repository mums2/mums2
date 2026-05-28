test_that("rarefy_ms returns the correct rowSum totals", {
  limit_cores()
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5),
                        number_of_threads = 2)
  results <- cluster_data(distances, dat,  0.3, "opticlust")
  community_object <- create_community_matrix_object(results)
  size <- 400
  resultant_matrix <- rarefy_ms(community_object, size, 10,
                                number_of_threads = 2)
  expect_true(all(rowSums(resultant_matrix) >= size))

})

test_that("rarefy_ms errors when given incorrect parameters", {
  limit_cores()
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5),
                        number_of_threads = 2)
  results <- cluster_data(distances, dat,  0.3, "opticlust")
  community_object <- create_community_matrix_object(results)
  expect_error(rarefy_ms("community_object", size, 10, 2),
               "Please ensure the community_object")
  expect_error(rarefy_ms(community_object, "size", 10, 2),
               "size must be numeric")
  expect_error(rarefy_ms(community_object, 1, "10", 2),
               "threshold must be numeric")
  expect_error(rarefy_ms(community_object, 1, 10, "a"),
               "number_of_threads must be numeric")
  expect_error(rarefy_ms(community_object, 1, 10, 1, "a"),
               "seed must be numeric")
})
