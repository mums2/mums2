#' @export
#' @title create community matrix
#' @description
#' Using your community_object, we are able to convert it into a community matrix for easier
#' usability of the object.
#' @param cluster_object the result of the `cluster_data()` function. 
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
#'   score_params = gnps_params(0.5), min_peaks = 0)
#' 
#' cluster_results <- cluster_data(distance_df = dist,
#'   ms2_match_data = matched_data, cutoff = 0.3, cluster_method = "opticlust")
#' 
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

#' @export
#' @title Convert Samples to Group Averages
#' @description
#' To account for users measuring there data in triplicates or other forms of measurement,
#' we have implemented a function that can tranform your matched data object to use
#' group averages instead of each sample individually. 
#' @param matched_data your mass data set object generated from `ms2_ms1_compare()`.
#' @param mpactr_object The object created from `import_all_data()`.
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
#' matched_data_avg <- convert_samples_to_group_averages(matched_data, filtered_data)
#' 
convert_samples_to_group_averages <- function(matched_data, mpactr_object) {
  trips <- t(get_triplicate_averages(mpactr_object, matched_data))
  meta_data <- get_meta_data(mpactr_object)
  injection_samples <- meta_data$Injection
  modifed_peak_table <- matched_data$ms1_data[,which(!(colnames(matched_data$ms1_data) %in% injection_samples)), with = FALSE]
  matched_data$ms1_data <- cbind(modifed_peak_table, trips)
  matched_data$samples <- unique(meta_data$Sample_Code)
  return(matched_data)
}

# Helper function for creating count tables
create_count_table <- function(ms2_match_data) {
  ms2_matches_compounds <- ms2_match_data$ms2_matches$ms1_compound_id
  peak_table <- ms2_match_data$ms1_data[ ,c("Compound", ms2_match_data$samples), with = FALSE]

  samples <- peak_table[which(peak_table$Compound %in% ms2_matches_compounds), ]
  return(data.frame(Representative_Sequence = samples$Compound, 
                    total = rowSums(samples[,-1]), samples[,-1], check.names = FALSE))
}


# Helper function for getting a matrix that displays the average of the triplicates
get_triplicate_averages <- function(mpactr_data, matched_data) {
  compounds <- matched_data$ms2_matches$ms1_compound_id
  peak <- get_peak_table(mpactr_data)
  meta_data <- get_meta_data(mpactr_data)
  sample_codes <- unique(meta_data$Sample_Code)
  triplicate_averages <- matrix(0, nrow(peak), 0)
  rownames(triplicate_averages) <- peak$Compound
  for(sample in sample_codes) {
    means <- as.matrix(rowMeans(peak[, meta_data$Injection[which(meta_data$Sample_Code == sample)], with = FALSE]))
    triplicate_averages <- cbind(triplicate_averages, means)
  }
  colnames(triplicate_averages) <- sample_codes
  return(t(triplicate_averages))
}
