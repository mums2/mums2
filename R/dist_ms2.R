#' Calculate pairwise distance between MS/MS features.
#'
#' @description
#' `dist_ms2` calculates and stores all non-zero distance values above
#'  the user defined cutoff (default = 0.3).
#'
#' @details
#' This function takes a `mass_data` object as input and calculates distance
#'  between ms2 peaks. Currently, MS1 features without MS2 peaks returns
#'  no distance value. Distance can be calculated with method `"gnps"`
#'  or `"spectral_entropy"`. A sparse matrix is returned.
#'
#'
#' @param data the object generated from `ms2_ms1_compare()`.
#' @param cutoff The maximum distance value (`numeric`) to store a pairwise
#'  comparison. The default of .3 corresponds to a cosine score of .7,
#'  meaning pairs with a score of .7 or higher will be stored in the matrix.
#' @param precursor_thresh Precursor mz tolerance. MS2 scans with a
#'  difference in precursor mz less than or equal to this value will be scored.
#' @param score_params Parameters for scoring method to be applied.
#'  See [gnps_params()] and [spec_entropy_params()] for more details.
#' @param min_peaks the minimum number of peaks that need to be present before
#' you compare the ms2 spectra.
#' @examples
#' squid_data <- import_all_data(peak_table = mums2::mums2_example("squid_peak_table.csv"), 
#'                             meta_data = mums2::mums2_example("squid_meta_data.csv"), 
#'                              format = "None")
#'
#' squid_filter <- squid_data |>
#'    filter_peak_table(filter_mispicked_ions_parameters()) |>
#'    filter_peak_table(filter_cv_parameters(cv_threshold = 0.2, cv_param = "mean")) |>
#'    filter_peak_table(filter_group_parameters(group_threshold = 0.1, "Blanks")) |>
#'    filter_peak_table(filter_insource_ions_parameters())
#'
#'
#' matched_data <- ms2_ms1_compare(mums2_example("12152023_Coculture_with_new_JC1.gnps.mgf"),
#'  squid_filter, 2, 6)
#' 
#' dist_gnps <- dist_ms2(data = matched_data, cutoff = 0.3, precursor_thresh = 2,
#'  score_params = gnps_params(0.5), min_peaks = 0)
#' 
#' dist_entropy <- dist_ms2(data = matched_data, cutoff = 0.3, precursor_thresh = 2,
#'  score_params = spec_entropy_params(), min_peaks = 0)
#'
#' @return A sparse matrix of class `"data.frame"`
#' @export
dist_ms2 <- function(data, cutoff, precursor_thresh, score_params, min_peaks = 6) {
  UseMethod("dist_ms2", data)
}

#' @method dist_ms2 mass_data
#' @export
dist_ms2.mass_data <- function(data, cutoff, precursor_thresh, score_params, min_peaks = 6) {
  data_list <- list("pmz" = data$ms2_matches$mz,
                    "id" = data$ms2_matches$ms1_compound_id,
                    "spectra" = data$peak_data)

  dist <- distMS2(data_list, score_params, precursor_thresh, cutoff, min_peaks)

  return(dist)
}
