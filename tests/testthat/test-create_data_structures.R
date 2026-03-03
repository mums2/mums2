test_that("test that we can create a community matrix", {
  dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5), min_peaks = 0)
  results <- cluster_data(distances, dat,  0.3, "opticlust")
  community_object <- create_community_matrix_object(results)
  mat <- create_community_matrix(results)
  expect_true("matrix" %in% class(mat))
  expect_true(all(mat == get_community_matrix(community_object)))
  expect_true(nrow(mat) == 21)
  expect_true(ncol(mat) == 339)
})

test_that("test that create a community matrix errors
          when given wrong inputs", {
            dat <- readRDS(test_path("exttestdata", "matched_data.RDS"))
            distances <- dist_ms2(dat, 0.3, 2, modified_cosine_params(0.5),
                                  min_peaks = 0)
            results <- cluster_data(distances, dat,  0.3, "opticlust")
            community_object <- create_community_matrix_object(results)
            expect_error(create_community_matrix(community_object))
          })


test_that("test that we create a proper count table", {
  ms2_match_data <- readRDS(test_path("exttestdata", "matched_data.RDS"))
  count_table <- create_count_table(ms2_match_data)
  ms2_matches_compounds <- ms2_match_data$ms2_matches$ms1_compound_id
  peak_table <- ms2_match_data$ms1_data[, -c(2, 3, 4)]
  peak_table$cor <- NULL
  samples <- peak_table[which(peak_table$Compound %in% ms2_matches_compounds), ]
  row_sums <- rowSums(samples[, -1])
  expect_true(all(names(count_table)[3:length(names(count_table))] %in%
                    names(peak_table[, -1])))
  expect_true(length(names(count_table)) == 23)
  expect_true(all(count_table$total == row_sums))
})


test_that("We can convert a mass_data object to an averaged mass data object", {
  data <-
    import_all_data(peak_table = test_path("exttestdata", "peak_table.csv"),
                    meta_data = test_path("exttestdata", "meta_data.csv"),
                    format = "Progenesis")

  mgf_files <- test_path("exttestdata",
                         "12152023_Coculture_with_new_JC1.gnps.mgf")
  ms2_data <- ms2_ms1_compare(mgf_files, data, 2, 6)
  ms2_avg_data <- convert_to_group_averages(ms2_data, data)
  meta_data <- get_meta_data(data)
  expect_true(all(meta_data$Sample_Code %in% ms2_avg_data$samples))
  expect_true(all(meta_data$Sample_Code %in% colnames(ms2_avg_data$ms1_data)))
})

test_that("convert to group averages fail if given wrong parameters", {
  data <-
    import_all_data(peak_table = test_path("exttestdata", "peak_table.csv"),
                    meta_data = test_path("exttestdata", "meta_data.csv"),
                    format = "Progenesis")

  mgf_files <- test_path("exttestdata",
                         "12152023_Coculture_with_new_JC1.gnps.mgf")
  ms2_data <- ms2_ms1_compare(mgf_files, data, 2, 6)
  expect_error(convert_to_group_averages("ms2_data", data),
               "The mass_data object")
  
  expect_error(convert_to_group_averages(ms2_data, "data"),
               "The mpactr object")
})



test_that("get_triplicate_averages returns a dataframe with all
          the triplicate averages", {
            data <-
              import_all_data(peak_table = test_path("exttestdata",
                                                     "peak_table.csv"),
                              meta_data = test_path("exttestdata",
                                                    "meta_data.csv"),
                              format = "Progenesis")

            mgf_files <- test_path("exttestdata",
                                   "12152023_Coculture_with_new_JC1.gnps.mgf")
            ms2_data <- ms2_ms1_compare(mgf_files, data, 2, 6)
            triplicate_data <-
              as.data.frame(get_triplicate_averages(data, ms2_data))
            meta_data <- get_meta_data(data)
            expect_true(all(rownames(triplicate_data)
                            == unique(meta_data$Sample_Code)))
            expect_true(ncol(triplicate_data) == nrow(get_peak_table(data)))
          })

test_that("generate_a_combined_table returns a data.frame with proper data", {
  data <-
    import_all_data(peak_table = test_path("exttestdata",
                                           "peak_table.csv"),
                    meta_data = test_path("exttestdata",
                                          "meta_data.csv"),
                    format = "Progenesis")

  mgf_files <- test_path("exttestdata",
                         "12152023_Coculture_with_new_JC1.gnps.mgf")
  ms2_data <- ms2_ms1_compare(mgf_files, data, 2, 6)
  matched_data_only <- generate_a_combined_table(matched_data = ms2_data)
  expect_true(ncol(matched_data_only) == 22)
  expect_true(all(ms2_data$samples %in% colnames(matched_data_only)))

  distances <- dist_ms2(ms2_data, 0.3, 2, modified_cosine_params(0.5),
                        min_peaks = 0)
  cluster_results <- cluster_data(distances, ms2_data,  0.3, "opticlust")
  matched_data_and_cluster_data <-
    generate_a_combined_table(matched_data = ms2_data,
                              cluster_data = cluster_results)
  expect_true(ncol(matched_data_and_cluster_data) == 23)
  expect_true(all(ms2_data$samples %in%
                    colnames(matched_data_and_cluster_data)))
  expect_true(all(which(matched_data_and_cluster_data$omus != "")
                  == which(matched_data_and_cluster_data$ms2_id != "")))

  psu_msmls <- read_msp(test_path("exttestdata", "database_data/PSU-MSMLS.msp"))
  annotations <- annotate_ms2(ms2_data, psu_msmls,
                              modified_cosine_params(0.5), 2000, 0, 0,
                              min_peaks = 0)

  all_data <- generate_a_combined_table(matched_data = ms2_data,
                                        annotations = annotations,
                                        cluster_data = cluster_results)

  expect_true(ncol(all_data) == 24)
  expect_true(all(ms2_data$samples %in% colnames(all_data)))
  expect_true(length(which(all_data$annotations != "")) > 0)


  data <- change_rt_to_seconds_or_minute(data, "seconds")
  ms2_data <- ms2_ms1_compare(mgf_files, data, 2, 100)
  with_rt_in_seconds_column <-
    generate_a_combined_table(matched_data = ms2_data)
  expect_true("RTINSECONDS" %in% colnames(with_rt_in_seconds_column))

  with_rt_in_seconds_column <-
    generate_a_combined_table(matched_data = ms2_data,
                              annotations, cluster_results)
  expect_true("RTINSECONDS" %in% colnames(with_rt_in_seconds_column))

  data <- change_rt_to_seconds_or_minute(data, "minutes")
  ms2_data <- ms2_ms1_compare(mgf_files, data, 2, 100)
  with_rt_in_minutes_column <-
    generate_a_combined_table(matched_data = ms2_data)
  expect_true("RTINMINUTES" %in% colnames(with_rt_in_minutes_column))
})


test_that("generate_a_combined_table will fail if sent the wrong parameters", {

  data <-
    import_all_data(peak_table = test_path("exttestdata",
                                           "peak_table.csv"),
                    meta_data = test_path("exttestdata",
                                          "meta_data.csv"),
                    format = "Progenesis")

  mgf_files <- test_path("exttestdata",
                         "12152023_Coculture_with_new_JC1.gnps.mgf")
  ms2_data <- ms2_ms1_compare(mgf_files, data, 2, 6)
  expect_error(generate_a_combined_table(matched_data = mgf_files),
               "matched_data must")
  expect_error(generate_a_combined_table(matched_data = ms2_data,
                                         annotations = mgf_files),
               "annotations must")
  expect_error(generate_a_combined_table(matched_data = ms2_data,
                                         cluster_data = mgf_files),
               "cluster_data must")

  psu_msmls <- read_msp(test_path("exttestdata", "database_data/PSU-MSMLS.msp"))
  annotations <- annotate_ms2(ms2_data, psu_msmls,
                              modified_cosine_params(0.5),
                              2000, 0, 0, min_peaks = 0)

  expect_error(generate_a_combined_table(matched_data = ms2_data,
                                         annotations = annotations[, 1:5]),
               "annotations must contain a column named 'name'")
  expect_error(generate_a_combined_table(matched_data = ms2_data,
                                         annotations = annotations
                                         [, c("name", "query_ms2_id")]),
               "annotations must contain a column named 'query_ms1_id'")
})
