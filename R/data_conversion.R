#' @export
#' @title Convert distance data frame to a distance object.
#' @description Converts the object generated from `dist_shared()` to a 
#' `dist` object.
#' @param distance_data_frame the object generated from the `dist_shared()` function.
#' @examples
#' squid_data <- import_all_data(peak_table = mums2::example("squid_peak_table.csv"), 
#'                             meta_data = mums2::example("squid_meta_data.csv"), 
#'                              format = "None")
#'
#' squid_filter <- squid_data |>
#'    filter_peak_table(filter_mispicked_ions_parameters()) |>
#'    filter_peak_table(filter_cv_parameters(cv_threshold = 0.2, cv_param = "mean")) |>
#'    filter_peak_table(filter_group_parameters(group_threshold = 0.1, "Blanks")) |>
#'    filter_peak_table(filter_insource_ions_parameters())
#'
#'
#' matched_data <- ms2_ms1_compare(example("12152023_Coculture_with_new_JC1.gnps.mgf"),
#'  squid_filter, 2, 6)
#' 
#' dist <- dist_ms2(data = matched_data, cutoff = 0.3, precursor_thresh = 2,
#'  score_params = gnps_params(0.5), min_peaks = 0)
#' 
#' cluster_results <- cluster_data(distance_df = dist,
#'  ms2_match_data = matched_data, cutoff = 0.3, cluster_method = "opticlust")
#' 
#' community_object <- create_community_matrix_object(cluster_results)
#' 
#' avg <- dist_shared(community_object = community_object, size = 4000,
#'  threshold = 100, diversity_index = "bray", iterations = 1)
#' 
#' community_object_to_distance_object(avg)
#' @return a `dist` object.
community_object_to_distance_object <- function(distance_data_frame) {
  dist_object <- matrix(distance_data_frame$diversity,
         sqrt(nrow(distance_data_frame)), sqrt(nrow(distance_data_frame)))
  samples <- unique(c(distance_data_frame$firstSample, distance_data_frame$otherSample))
  rownames(dist_object) <- samples
  colnames(dist_object) <- samples
  return(as.dist(dist_object))
}