test_that("test that we can create a community matrix", {
  count_table <- test_path("exttestdata", "final.count_table")
  distances <- test_path("exttestdata", "final.dist")
  final_count <- read_count(count_table)
  final_dist <- read_dist(distances, final_count, 0.03, F)
  final_cluster <- cluster(final_dist, 0.03, "opticlust")
  mat <- create_community_matrix(final_cluster)
  communiy_object <- create_community_matrix_object(final_cluster)
  expect_true("matrix" %in% class(mat))
  expect_true(all(mat == get_community_matrix(communiy_object)))
  expect_true(nrow(mat) == 19)
  expect_true(ncol(mat) == 484)
})

test_that("test that we create a proper count table", {
  data <- import_all_data(peak_table = test_path("exttestdata", "peak_table.csv"), 
  meta_data = test_path("exttestdata", "meta_data.csv"), 
  format = "Progenesis")

  data_filtered <- data |>
    filter_peak_table(filter_mispicked_ions_parameters())

  mass_data_set <- convert_mpactr_object_to_mass_data_set(data_filtered)
  count_table <- create_count_table(mass_data_set)
  row_sums <- rowSums(mass_data_set@expression_data)
  expect_true(all(names(count_table)[3:length(names(count_table))] %in% 
                  names(mass_data_set@expression_data)))
  expect_true(length(names(count_table)) == 20)
  expect_true(all(count_table$total == row_sums))

})
