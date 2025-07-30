test_that("test that we can create a community matrix", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5), min_peaks = 0)
  results <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(results)
  mat <- create_community_matrix(results)
  expect_true("matrix" %in% class(mat))
  expect_true(all(mat == get_community_matrix(communiy_object)))
  expect_true(nrow(mat) == 21)
  expect_true(ncol(mat) == 339)
})

test_that("test that create a community matrix errors when given wrong inputs", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, gnps_params(0.5), min_peaks = 0)
  results <- cluster_data(distances, dat,  0.3, "opticlust")
  communiy_object <- create_community_matrix_object(results)
  expect_error(create_community_matrix(communiy_object))
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


test_that("We can convert a mass_data object to an averaged mass data object",
{
  data <- import_all_data(peak_table = test_path("exttestdata", "peak_table.csv"), 
                          meta_data = test_path("exttestdata", "meta_data.csv"), 
                          format = "Progenesis")
                          
  mgf_files <- test_path("exttestdata", "12152023_Coculture_with_new_JC1.gnps.mgf")
  ms2_data <- ms2_ms1_compare(mgf_files, data, 2, 6)
  ms2_avg_data <- convert_samples_to_group_averages(ms2_data, data)
  meta_data <- get_meta_data(data)
  expect_true(all(meta_data$Sample_Code %in% ms2_avg_data$samples))
  expect_true(all(meta_data$Sample_Code %in% colnames(ms2_avg_data$ms1_data)))
})


test_that("get_triplicate_averages returns a dataframe with all the triplicate averages",
{
  data <- import_all_data(peak_table = test_path("exttestdata", "peak_table.csv"), 
                          meta_data = test_path("exttestdata", "meta_data.csv"), 
                          format = "Progenesis")
                          
  mgf_files <- test_path("exttestdata", "12152023_Coculture_with_new_JC1.gnps.mgf")
  ms2_data <- ms2_ms1_compare(mgf_files, data, 2, 6)
  triplicate_data <- as.data.frame(get_triplicate_averages(data, ms2_data))
  meta_data <- get_meta_data(data)
  expect_true(all(rownames(triplicate_data) == unique(meta_data$Sample_Code)))
  expect_true(ncol(triplicate_data) == nrow(get_peak_table(data)))
})
