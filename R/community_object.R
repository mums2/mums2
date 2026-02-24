#' @export
#' @title Create Community Matrix Object.
#' @description
#' Using the data generated from clustering or adding ms2 data to your object,
#' we are able to create a community matrix object. The community matrix object
#' stores the same data a community matrix but within a cpp object. We use this object
#' to conduct analysis more efficiently. 
#' @param data the result of the `cluster_data()` function,
#'  or just a mass_data object created from `ms2_ms1_compare()`.
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
#' community_with_cluster <- create_community_matrix_object(cluster_results)
#'
#' community_object_mass_data <- create_community_matrix_object(matched_data)
#'
#' @return a external pointer to an Rcpp object.
create_community_matrix_object <- function(data) {
  if(!inherits(data, "mass_data") &&
     !inherits(data, "mothur_cluster")) {
    stop(paste0("data has to be created using the",
                " `ms2_ms1_compare()` or `cluster_data()` functions"))
  }
  UseMethod("create_community_matrix_object", data)
}

#' @export
#' @rdname create_community_matrix_object
create_community_matrix_object.mass_data <- function(data) {
  samples <- data$samples
  ms2_matches <-  data$ms2_matches$ms1_compound_id
  filtered_data <- data$ms1_data[which(data$ms1_data$Compound %in% ms2_matches),
  ][, samples, with = FALSE]
  matrix <- as.matrix(t(filtered_data))
  rownames(matrix) <- samples
  colnames(matrix) <- ms2_matches
  community_matrix <- CreateCommunityMatrix(matrix)
  class(community_matrix) <- c(class(community_matrix), "community_object")
  community_matrix
}

#' @export
#' @rdname create_community_matrix_object
create_community_matrix_object.mothur_cluster <- function(data) {
  df <- data$abundance
  samples <- unique(df$samples)
  combined_df <- data.frame(abund = df[which(df$samples == samples[[1]]),
                            ]$abundance)

  for (i in 2:length(samples)) {
    combined_df <- cbind(combined_df,
                         data.frame(abund = df[which(df$samples == samples[[i]])
                                               , ]$abundance))
  }

  combined_df <- t(as.matrix(combined_df))
  rownames(combined_df) <- samples
  colnames(combined_df) <- data$cluster$omu
  obj <- CreateCommunityMatrix(combined_df)
  class(obj) <- c(class(obj), "community_object")
  obj
}

#' @export
#' @title Get Community Matrix
#' @description
#' Returns the community `matrix` or the data
#'  that you used to create the object.
#' @param community_object the object created from
#'  the `create_community_object()` function.
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
#' community_with_cluster <- create_community_matrix_object(cluster_results)
#'
#' community_object_mass_data <- create_community_matrix_object(matched_data)
#' @return returns `matrix`, based on the community object.
get_community_matrix <- function(community_object) {
  if (!("community_object" %in% class(community_object))) {
    stop("Please ensure the community_object is created
         from the `create_community_object` function.")
  }
  GetCommunityMatrix(community_object)
}

#' @export
#' @title Print Community Object
#' @description
#' S3 function for print the community object
#' @param x the object created from the `create_community_object()` function.
#' @param ... other parameters that are included in the print function.
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
#' community_object <- create_community_matrix_object(cluster_results)
#' print(community_object)
#'
print.community_object <- function(x, ...) {
  print(get_community_matrix(x), ...)
}
