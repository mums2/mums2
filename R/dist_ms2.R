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
#' @param precursor_threshold Precursor mz tolerance. MS2 scans with a
#'  difference in precursor mz less than or equal to this value will be scored.
#'  Disable this by setting this value to -1 or less.
#' @param score_params Parameters for scoring method to be applied.
#'  See [modified_cosine_params()] and [spec_entropy_params()] for more details.
#' @param min_peaks the minimum number of peaks that need to be present before
#' you compare the ms2 spectra.
#' @param number_of_threads the number of
#' threads you want to use for this calculation.
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
#'
#'
#' matched_data <- ms2_ms1_compare(mums2_example("botryllus_v2.gnps.mgf"),
#'  filtered_data, 10, 6)
#'
#' dist_gnps <- dist_ms2(data = matched_data,
#'  cutoff = 0.3, precursor_threshold = 2,
#'  score_params = modified_cosine_params(0.5), min_peaks = 0)
#'
#' dist_entropy <- dist_ms2(data = matched_data,
#'  cutoff = 0.3, precursor_threshold = 2,
#'  score_params = spec_entropy_params(), min_peaks = 0)
#'
#' @return A sparse matrix of class `"data.frame"`
#' @export
dist_ms2 <- function(data, cutoff, precursor_threshold, score_params,
                     min_peaks = 6,
                     number_of_threads = detectCores()) {
  if (!inherits(data, "mass_data")) {
    stop(paste0("The mass_data object must be created using the",
                " `ms2_ms1_compare() function`"))
  }
  if (nrow(data$ms2_matches) <= 0) {
    stop("Cannot calculate distances, there are no matched ms2.")
  }
  if (!inherits(score_params, "parameters")) {
    stop(paste0("score_params must be created using the",
                " modified_cosine_params() or spec_entropy_params() function"))
  }

  if (!is.numeric(cutoff)) {
    stop("cutoff must be numeric")
  }

  if (!is.numeric(precursor_threshold)) {
    stop(paste0("precursor_threshold must be a numeric"))
  }

  if (!is.numeric(min_peaks)) {
    stop("min_peaks must be a numeric")
  }

  if (!is.numeric(number_of_threads)) {
    stop("number_of_threads must be a numeric")
  }


  data_list <- list("pmz" = data$ms2_matches$mz,
                    "id" = data$ms2_matches$ms1_compound_id,
                    "spectra" = data$peak_data)

  dist <- distMS2(data_list, score_params, precursor_threshold,
                  cutoff, min_peaks, number_of_threads)

  class(dist) <- c(class(dist), "mass_data_dist")
  return(dist)
}
