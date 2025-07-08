#' @export
#' @title create community matrix
#' @description
#' Using your community_object, we are able to convert it into a community matrix for easier
#' usability of the object.
#' @param cluster_object the result of the `cluster_data()` function. 
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
#'   squid_filter, 2, 6)
#' dist <- dist_ms2(data = matched_data, cutoff = 0.3, precursor_thresh = 2,
#'   score_params = gnps_params(0.5), min_peaks = 0)
#' cluster_results <- cluster_data(distance_df = dist,
#'   ms2_match_data = matched_data, cutoff = 0.3, cluster_method = "opticlust")
#' community_matrix <- create_community_matrix_object(cluster_results)
#' @return a `data.frame` object of your community_object.
create_community_matrix <- function(cluster_object) {
  df <- get_abundance(cluster_object)
  samples <- unique(df$samples)
  combined_df <- data.frame(abund = df[which(df$samples == samples[[1]]), ]$abundance)

  for(i in 2:length(samples)) {
    combined_df <- cbind(combined_df, data.frame(abund = df[which(df$samples == samples[[i]]), ]$abundance))
  }

  combined_df <- t(as.matrix(combined_df))
  rownames(combined_df) <- samples
  return(combined_df)
}

# Helper function for creating count tables
create_count_table <- function(ms2_match_data) {
  ms2_matches_compounds <- ms2_match_data$ms2_matches$ms1_compound_id
  peak_table <- ms2_match_data$ms1_data[ ,c("Compound", ms2_match_data$samples), with = FALSE]

  samples <- peak_table[which(peak_table$Compound %in% ms2_matches_compounds), ]
  return(data.frame(Representative_Sequence = samples$Compound, 
                    total = rowSums(samples[,-1]), samples[,-1], check.names = FALSE))
  }