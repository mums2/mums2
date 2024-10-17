#' Calculate pairwise distance between MS/MS features.
#'
#' @description
#' `dist_ms2` calculates and stores all non-zero distance values above
#'  the user defined cutoff (default = 0.3).
#'
#' @details
#' This function takes a mass_dataset as input and calculates distance
#'  between ms2 peaks. Currently, MS1 features without MS2 peaks returns
#'  no distance value. Distance can be calculated with method `"gnps"`
#'  or `"spectral_entropy"`. A sparse matrix is returned.
#'
#'
#' @param data An object containing MA2 data.
#' @param cutoff The maximum distance value (`numeric`) to store a pairwise
#'  comparison. The default of .3 corresponds to a cosine score of .7,
#'  meaning pairs with a score of .7 or higher will be stored in the matrix.
#' @param precursor_thresh Precursor mz tolerance. MS2 scans with a
#'  difference in precursor mz less than or equal to this value will be scored.
#' @param score_params Parameters for scoring method to be applied.
#'  See [gnps_params()] and [spec_entropy_params()] for more details.
#'
#' @return A sparse matrix of class `"data.frame"`
#' @export
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#'
#' dat <- tidyMassDemo
#' first_5 <- c("M71T776_POS", "M72T54_POS", "M81T51_POS",
#'    "M83T51_1_POS", "M83T51_2_POS")
#'
#' dat_sub <- dat %>%
#'    massdataset::activate_mass_dataset("variable_info") %>%
#'    massdataset::filter(variable_id %in% first_5)
#'
#' dat_sub_dist <- dist_ms2(dat_sub, 0.3, 2, gnps_params(0.5))}
dist_ms2 <- function(data, cutoff, precursor_thresh, score_params) {
  UseMethod("dist_ms2", data)
}

#' @method dist_ms2 mass_dataset
#' @export
dist_ms2.mass_dataset <- function(data, cutoff, precursor_thresh, score_params) {
  # ms1_data <- data@variable_info
  # n <- nrow(ms1_data)
  # ms2_matches <- data@ms2_data[[1]]@variable_id
  n <- length(data@ms2_data[[1]]@variable_id)

  dist <- MatrixExtra::emptySparse(nrow = n, ncol = n, format = "T", dtype = "d")

  if (any(duplicated(data@ms2_data[[1]]@variable_id))) {
    cli::cli_abort(c("One or more ms1 features are represented multiple times in ms2 scans."))
  }

  for (i in seq_along(1:(n-1))) {
    # if there is no ms2 scan for i
    # if (isFALSE(does_ms1_have_ms2(data$variable_id[i], 
    #                               ms2_matches))) {
    #   next
    # }
    # for(int i =0; )
    for (j in seq(i + 1, n)[seq(i + 1, n) <= n]) {
      val1 <- data@ms2_data[[1]]@ms2_mz[[i]]
      if (data@ms2_data[[1]]@ms2_mz[[i]] == data@ms2_data[[1]]@ms2_mz[[j]] &
          nrow(data@ms2_data[[1]]@ms2_spectra[[i]] ==
          nrow(data@ms2_data[[1]]@ms2_spectra[[j]]))) {
        dist[i, j] <- .Machine$double.xmin
        next
      }
      
      if (abs(data@ms2_data[[1]]@ms2_mz[[i]] - data@ms2_data[[1]]@ms2_mz[[j]]) > precursor_thresh) {
        # dist[i, j] <- .Machine$double.xmin
        next
      }

      score <- score_ms2(get_peaks_data(data, i), get_peaks_data(data, j),
                         score_params)
      
      if ((1 - score) <= cutoff) {
        dist[i, j] <- 1 - score
      }
    }

      
  }
  
  return(dist)
}

dist_ms2_cpp <- function(data, cutoff, precursor_thresh, score_params) {
  # spectraDataList, parameters, precursor_thresh, cutoff
  data_list <- list("pmz" = data@ms2_data[[1]]@ms2_mz,
                    "id" = data@ms2_data[[1]]@variable_id,
                    "spectra" = lapply(data@ms2_data[[1]]@ms2_spectra,
                                       as.data.frame))
  
  dist <- distMS2(data_list, score_params, precursor_thresh, cutoff)

  return(dist)
}

# is_same_scan <- function(spectra_1, spectra_2) {
#   if (nrow(spectra_1) != nrow(sepectra_2)) {
#     return(FALSE)
#   }

#   return(TRUE)
# }