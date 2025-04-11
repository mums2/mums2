test_that("test that we can create a community matrix", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5))
  results <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(results)
  mat <- get_community_matrix(communiy_object)
  expect_true("matrix" %in% class(mat))
  expect_true(all(mat == get_community_matrix(communiy_object)))
  expect_true(nrow(mat) == 21)
  expect_true(ncol(mat) == 214)
})

test_that("test that we create a proper count table", {
  ms2_match_data <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  count_table <- create_count_table(ms2_match_data)
  ms2_matches_compounds <- ms2_match_data$ms2_matches$ms1_compound_id
  peak_table <- ms2_match_data$ms1_data[ ,-c(2, 3, 4)]
  peak_table$cor <- NULL
  samples <- peak_table[which(peak_table$Compound %in% ms2_matches_compounds), ]
  row_sums <- rowSums(samples[,-1])
  expect_true(all(names(count_table)[3:length(names(count_table))] %in% 
                  names(peak_table[, -1])))
  expect_true(length(names(count_table)) == 23)
  expect_true(all(count_table$total == row_sums))
})

