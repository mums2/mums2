test_that("test that we can create a community matrix", {
  mass_data_set <- readRDS(test_path("exttestdata", "test_mass_data_set.RDS"))
  distances <- dist_ms2(mass_data_set, 0.3, 2, gnps_params(0.5))
  results <- cluster_data(distances, mass_data_set, "opticlust")
  communiy_object <- create_community_matrix_object(results)
  mat <- get_community_matrix(communiy_object)
  expect_true("matrix" %in% class(mat))
  expect_true(all(mat == get_community_matrix(communiy_object)))
  expect_true(nrow(mat) == 18)
  expect_true(ncol(mat) == 378)
})

test_that("test that we create a proper count table", {
  mass_data_set <- readRDS(test_path("exttestdata", "test_mass_data_set.RDS"))
  count_table <- create_count_table(mass_data_set)
  ms2_matches <- mass_data_set@ms2_data[[1]]@variable_id
  row_sums <- rowSums(mass_data_set@expression_data[which(rownames(mass_data_set@expression_data) %in% ms2_matches), ])
  expect_true(all(names(count_table)[3:length(names(count_table))] %in% 
                  names(mass_data_set@expression_data)))
  expect_true(length(names(count_table)) == 20)
  expect_true(all(count_table$total == row_sums))
})
