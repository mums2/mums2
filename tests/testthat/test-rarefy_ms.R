test_that("rarefy_ms returns the correct rowSum totals", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5))
  results <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(results)
  size <- 400
  resultant_matrix <- rarefy_ms(communiy_object, size, 10)
  expect_true(all(rowSums(resultant_matrix) >= size))

})

test_that("rarefy_ms errors when given an incorrect object", {
  expect_error(rarefy_ms(c(), 100, 10))
})
