#' Rarefy MS1 Feature Table
#' 
#' @description
#' `rarefy_ms()` performs a single subsampling of MS1 features in sample.
#'  Feature intensities are subsampled to the supplied `size` and accounts
#'  for intesnity thresholds due to machine limits and background noise.
#'  Specifically, features whose abundance falls below the `threshold`
#'  after rarefying are removed. This allows for accurate represetation
#'  of samples at diffrent dilutions regardless of the desired
#'  submsampling `size`.
#' 
#' @param community_object A `community_object`
#' @param size The desired total sample intensity to subsample to.
#' @param threshold The individual feature threshold. Each subsampled feature
#'  must be >= this value to be retained.
#' @return A `external_pointer` that references a community matrix of rarefied feature intensities.
#' @export
#' 
#' @examples
#' squid_data <- import_all_data(peak_table = mums2::mums2_example("squid_peak_table.csv"), 
#'                             meta_data = mums2::mums2_example("squid_meta_data.csv"), 
#'                              format = "None")
#'
#' squid_filter <- squid_data |>
#'    filter_peak_table(filter_mispicked_ions_parameters()) |>
#'    filter_peak_table(filter_cv_parameters(cv_threshold = 0.2)) |>
#'    filter_peak_table(filter_group_parameters(group_threshold = 0.1, "Blanks")) |>
#'    filter_peak_table(filter_insource_ions_parameters())
#'
#'
#' matched_data <- ms2_ms1_compare(mums2_example("12152023_Coculture_with_new_JC1.gnps.mgf"),
#'  squid_filter, 2, 6)
#' 
#' dist <- dist_ms2(data = matched_data, cutoff = 0.3, precursor_thresh = 2,
#'  score_params = gnps_params(0.5), min_peaks = 0)
#' 
#' cluster_results <- cluster_data(distance_df = dist,
#'  ms2_match_data = matched_data, cutoff = 0.3, cluster_method = "opticlust")
#' 
#' community_object <- create_community_matrix_object(cluster_results)
#' rarefy_ms(community_object, 4000, 100)
#' 
rarefy_ms <- function(community_object, size, threshold) {
  if(!("community_object" %in% class(community_object))) {
    stop("Please ensure the community_object is created from the `create_community_object` function.")
  }
  return(RarefactionCalculation(community_object, size, threshold))
}
