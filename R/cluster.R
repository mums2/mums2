#' @export
#' @title Cluster Features
#' @description
#' `cluster_data()` allows users to cluster features inside
#' the mass data object. This is done by creating a sparse matrix
#' using the `distMs2()` function and inputting that inside the
#' clutur package. This allows us to easily cluster features
#' that contain an ms2 spectra.
#' @param distance_df a distance df that was generated
#' from the `distMs2()` function.
#' @param ms2_match_data your mass data set object generated
#'  from `ms2_ms1_compare()`.
#' @param cutoff the cutoff value you wish to cluster to.
#' @param cluster_method the clustering algorithm you wish to use.
#'  The options are: furthest, nearest, weighted, average, and opticlust.
#' @examples
#' data <-
#'    import_all_data(peak_table =
#'                    mums2::mums2_example("botryllus_pt_small.csv"),
#'                    meta_data =
#'                    mums2::mums2_example("meta_data_boryillus.csv"),
#'                    format = "None")
#'
#' filtered_data <- data |>
#'    filter_peak_table(filter_mispicked_ions_params()) |>
#'    filter_peak_table(filter_cv_params(cv_threshold = 0.2)) |>
#'    filter_peak_table(filter_group_params(group_threshold = 0.1,
#'                                              "Blanks")) |>
#'    filter_peak_table(filter_insource_ions_params())
#'
#' change_rt_to_seconds_or_minute(filtered_data, "minutes")
#'
#' matched_data <- ms2_ms1_compare(mums2_example("botryllus_v2.gnps.mgf"),
#'  filtered_data, 10, 6)
#'
#' dist <- dist_ms2(data = matched_data, cutoff = 0.3, precursor_thresh = 2,
#'  score_params = modified_cosine_params(0.5), min_peaks = 0)
#'
#' cluster_results <- cluster_data(distance_df = dist,
#'  ms2_match_data = matched_data, cutoff = 0.3, cluster_method = "opticlust")
#'
#' @return a shared `data.frame` (or a `mothur_cluster` object) displaying all
#' the clustered and abundance data.
cluster_data <- function(distance_df, ms2_match_data,
                         cutoff = 0.3, cluster_method = "opticlust") {

  if (!inherits(distance_df, "mass_data_dist")) {
    stop("distance_df must be an object created from the `dist_ms2()` function")
  }

  if (nrow(distance_df) <= 0) {
    stop("distance_df must have more than 0 rows")
  }

  if (!inherits(ms2_match_data, "mass_data")) {
    stop(paste0("The mass_data object must be created using the",
                " `ms2_ms1_compare()`"))
  }

  if (!is.numeric(cutoff)) {
    stop("cutoff should be a numeric value")
  }
  sparse_matrix <- create_sparse_matrix(distance_df$i,
                                        distance_df$j, distance_df$dist)
  # Create Count Table
  count_table <- create_count_table(ms2_match_data)

  # Create Distance Object
  dist <- read_dist(sparse_matrix, count_table, cutoff, FALSE)

  # Cluster Data
  return(cluster(dist, cutoff, cluster_method, bin_column_name_to = "omu"))
}
