#' @title Distance Shared
#' @export
#' @description dissimilarity via beta diversity
#' @param community_object the object created from
#' the `create_community_object()` function.
#' @param size the size you wish to rarefy your diversity matrix to.
#' @param threshold the threshold you want your species to
#' reach before it is included in the rarefaction sum.
#' @param diversity_index the diversity index you wish
#' to calculate diversity. You can choose from:
#' bray, jaccard, soren, hamming, morista, and thetayc.
#' @param subsample if true, we will rarefy the data before we run
#' the diversity calculations. Default is TRUE.
#' @param number_of_threads the amount of threads you want the
#' calculation to use.
#' @param iterations the amount of times you wish to run your calculation.
#' @param seed the RNG (random number generator) seed you would like to use.
#' @examples
#' data <- import_all_data(peak_table =
#'                         mums2::mums2_example("full_mix_peak_table.csv"),
#'                         meta_data =
#'                         mums2::mums2_example("full_mix_meta_data.csv"),
#'                         format = "Metaboscape")
#'
#' filtered_data <- data |>
#'    filter_peak_table(filter_mispicked_ions_parameters()) |>
#'    filter_peak_table(filter_cv_parameters(cv_threshold = 0.2)) |>
#'    filter_peak_table(filter_group_parameters(group_threshold = 0.1,
#'                                              "Blanks")) |>
#'    filter_peak_table(filter_insource_ions_parameters())
#'
#'  change_rt_to_seconds_or_minutes(filtered_data, "minutes")
#'
#' matched_data <- ms2_ms1_compare(mums2_example("full_mix_ms2.mgf"),
#'  filtered_data, 2, 6)
#' dist <- dist_ms2(data = matched_data, cutoff = 0.3, precursor_thresh = 2,
#'  score_params = gnps_params(0.5), min_peaks = 0)
#'
#' cluster_results <- cluster_data(distance_df = dist,
#'  ms2_match_data = matched_data,
#'  cutoff = 0.3, cluster_method = "opticlust")
#'
#' community_object <- create_community_matrix_object(cluster_results)
#'
#' dist_shared(community_object, 4000, 100, "bray", TRUE, 1)
#' @return a `data.frame` object that shows the
#' dissimilarity between all samples.
dist_shared <- function(community_object, size, threshold,
                        diversity_index = "bray", subsample = TRUE,
                        number_of_threads = detectCores(),
                        iterations = 100, seed = 123) {
  diversity_index_list <- c("bray", "jaccard", "soren",
                            "hamming", "morisita", "thetayc")
  if (!("community_object" %in% class(community_object))) {
    stop("Please ensure the community_object is created from the
         `create_community_object` function.")
  }
  if (!(diversity_index %in% diversity_index_list)) {
    stop(paste0("Please ensure your diversity
                index is one of the following values: ",
                paste(diversity_index_list, collapse = ", "))
    )
  }
  result <- FasterAvgDist(community_object, diversity_index,
                          size, threshold, subsample,
                          number_of_threads, iterations, seed)
  result[which(is.nan(result))] <- 0
  return(as.dist(result))
}


#' @title Alpha Diversity Summary
#' @export
#' @description alpha diversity
#' @param community_object the object created from
#' the `create_community_object()` function.
#' @param size the size you wish to rarefy your diversity matrix to.
#' @param threshold the threshold you want your
#' species to reach before it is included in the rarefaction sum.
#' @param diversity_index the diversity index you wish to calculate diversity,
#' the two options are shannon or simpson.
#' @param subsample if true, we will rarefy the data before we
#' run the diversity calculations. Default is TRUE.
#' @param number_of_threads the amount of threads
#' you want the calculation to use.
#' @param iterations the amount of times you wish to run your calculation.
#' @param seed the RNG (random number generator) seed you would like to use.
#' @examples
#' #' data <- import_all_data(peak_table =
#'                         mums2::mums2_example("full_mix_peak_table.csv"),
#'                         meta_data =
#'                         mums2::mums2_example("full_mix_meta_data.csv"),
#'                         format = "Metaboscape")
#'
#' filtered_data <- data |>
#'    filter_peak_table(filter_mispicked_ions_parameters()) |>
#'    filter_peak_table(filter_cv_parameters(cv_threshold = 0.2)) |>
#'    filter_peak_table(filter_group_parameters(group_threshold = 0.1,
#'                                              "Blanks")) |>
#'    filter_peak_table(filter_insource_ions_parameters())
#'
#'  change_rt_to_seconds_or_minutes(filtered_data, "minutes")
#'
#' matched_data <- ms2_ms1_compare(mums2_example("full_mix_ms2.mgf"),
#'  filtered_data, 2, 6)
#' dist <- dist_ms2(data = matched_data, cutoff = 0.3, precursor_thresh = 2,
#'  score_params = gnps_params(0.5), min_peaks = 0)
#'
#' cluster_results <- cluster_data(distance_df = dist,
#'  ms2_match_data = matched_data,
#'  cutoff = 0.3, cluster_method = "opticlust")
#'
#' community_object <- create_community_matrix_object(cluster_results)
#'
#' alpha_summary(community_object, 4000, 100, "shannon", TRUE, iterations = 1)
#' @return a `data.frame` object that shows the dissimilarity in samples.
alpha_summary <- function(community_object, size, threshold,
                          diversity_index = "shannon",
                          subsample = TRUE,
                          number_of_threads = detectCores(),
                          iterations = 1000, seed = 123) {
  diversity_index_list <- c("shannon", "simpson")
  if (!("community_object" %in% class(community_object))) {
    stop("Please ensure the community_object is created from the
         `create_community_object` function.")
  }

  if (!(diversity_index %in% diversity_index_list)) {
    stop(paste0("Please ensure your diversity index
                is one of the following values: ",
                paste(diversity_index_list, collapse = ", "))
    )
  }

  result <- FasterAvgDist(community_object, diversity_index, size, threshold,
                          subsample, number_of_threads, iterations, seed)
  result[which(is.nan(result))] <- 0
  return(result)
}
