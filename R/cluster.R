#' @export
#' @title Cluster
#' @description
#' Clusters the data together
#' @param distance_df a distance df that was generated from the `distMs2()` function.
#' @param ms2_match_data your mass data set object generated from `ms2_ms1_compare()`.
#' @param cutoff the cutoff value you wish to cluster to.
#' @param cluster_method a cluster method, there are five methods to choose from:
#' furthest, nearest, weighted, average, and opticlust.
#' @examples 
#' data <- import_all_data(peak_table = mums2::mums2_example("full_mix_peak_table.csv"), 
#'                             meta_data = mums2::mums2_example("full_mix_meta_data.csv"), 
#'                              format = "Metaboscape")
#'
#' filtered_data <- data |>
#'    filter_peak_table(filter_mispicked_ions_parameters()) |>
#'    filter_peak_table(filter_cv_parameters(cv_threshold = 0.2)) |>
#'    filter_peak_table(filter_group_parameters(group_threshold = 0.1, "Blanks")) |>
#'    filter_peak_table(filter_insource_ions_parameters())
#'
#' change_rt_to_seconds_or_minutes(filtered_data, "minutes")
#' 
#' matched_data <- ms2_ms1_compare(mums2_example("full_mix_ms2.mgf"),
#'  filtered_data, 2, 6)
#' 
#' dist <- dist_ms2(data = matched_data, cutoff = 0.3, precursor_thresh = 2,
#'  score_params = gnps_params(0.5), min_peaks = 0)
#' 
#' cluster_results <- cluster_data(distance_df = dist,
#'  ms2_match_data = matched_data, cutoff = 0.3, cluster_method = "opticlust")
#'
#' @return a shared `data.frame` displaying all the clustered and abundance data.
cluster_data <- function(distance_df, ms2_match_data, cutoff = 0.3, cluster_method = "opticlust") {

  sparse_matrix <- create_sparse_matrix(distance_df$i, distance_df$j, distance_df$dist)
  # Create Count Table 
  count_table <- create_count_table(ms2_match_data)
  
  # Create Distance Object
  dist <- read_dist(sparse_matrix, count_table, cutoff, FALSE)

  # Cluster Data
  return(cluster(dist, cutoff, cluster_method, bin_column_name_to = "omu"))
}