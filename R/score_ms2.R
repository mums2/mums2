#' Calculate similarity score between two ms2 spectra
#'
#' @description
#' This function is calculates similarity between two ms2 spectra.
#'  Currently two methods for similarity are supported: `"gnps"` and
#'  `"spectral_entropy"`.
#'
#' Similaritity between level 2 spectra are determined via
#' spectral scoring methods. Currently scoring methods `"gnps"` and
#' `"spectral_entropy"` are supported. The scoring method is specified by the
#' `score_params` argument. `score_params` is a list of parameters for the
#' chosen scoring method. Parameters for "gnps" and "spectral_entropy" can be
#' created with notur functions [gnps_params()] and [spec_entropy_params()],
#' respectively.
#'
#' @details
#' [gnps_params()] will calculate GNPS cosine score, based on code by Wang
#' et al. (2016).This method applies a square root normalization to peak
#' intensities, aligns peaks both with and without correction for mass
#' shifts, and calculates similarity.
#'
#' [spec_entropy_params()] will calculate entropy similarity using functions
#'  from the `msentropy` package (Li et al. 2021) to calculate. For more
#' information about parameters, see [msentropy::msentropy_similarity()].
#'
#' @param peaks_1 A `data.frame` with two columns: `mz` and `intensity`.
#' @param peaks_2 A `data.frame` with two columns: `mz` and `intensity`.
#' @param score_params Parameters for scoring method to be applied.
#'
#' @return Similarity score.
#' @export
#'
#' @references
#' Mingxun Wang, Jeremy J. Carver, Vanessa V. Phelan, Laura M. Sanchez,
#'  Neha Garg, Yao Peng, Don Duy Nguyen et al. "Sharing and community
#'  curation of mass spectrometry data with Global Natural Products Social
#'  Molecular Networking." Nature biotechnology 34, no. 8 (2016): 828.
#'  PMID: 27504778
#'
#' Li, Y., Kind, T., Folz, J. et al. Spectral entropy outperforms MS/MS
#'  dot product similarity for small-molecule compound identification,
#'  Nat Methods 18, 1524–1531 (2021). https://doi.org/10.1038/s41592-021-01331-z
#'
#' @examples
#' pmz_1 <- 136.0615
#' mz_1 <- c(53.0124, 53.0370, 55.0283, 56.0316, 65.0128, 65.0373, 66.0157,
#'           66.0337, 67.0282, 68.0323,  68.0467, 77.0118)
#' int_1 <- c(1000, 290, 6800, 200, 16000, 410, 1200, 190, 12000, 600,
#'                  320, 1800)
#' peaks_1 <- create_peaks_data(data.frame("mz" = mz_1,
#'                                          "intensity" = int_1),
#'                               pmz_1)
#'
#'
#' pmz_2 <- 136.062
#' mz_2 <- c(52.0220, 53.0349, 53.9952, 53.9984, 54.0015, 54.0025, 54.0046,
#'           54.0098, 54.0243, 54.0264, 54.8659, 55.0564, 55.9345, 55.9355,
#'           65.0155, 66.0233, 66.0279, 66.0382, 66.8263, 66.8298, 67.0319,
#'           67.0492, 67.0515, 67.0538, 67.0585, 68.0272, 71.9321, 72.9389,
#'           72.9413)
#' int_2 <- c(8, 276, 15, 6, 7, 753, 21, 15, 9, 25, 29, 11, 6, 22, 44, 33,
#'            5, 3959, 41, 24, 18, 9, 171, 91, 5, 19, 13, 1481, 13)
#' peaks_2 <- create_peaks_data(data.frame("mz" = mz_2,
#'                                          "intensity" = int_2),
#'                               pmz_2)
#'
#' score <- score_ms2(peaks_1, peaks_2, gnps_params(frag_tolerance = 0.5))
#' score <- score_ms2(peaks_1, peaks_2, spec_entropy_params())
#'
score_ms2 <- function(peaks_1, peaks_2, score_params) {
  UseMethod("score_ms2", peaks_1)
}


#' @method score_ms2 peaks_data
#' @export
score_ms2.peaks_data <- function(peaks_1, peaks_2, score_params) {
  score <- ScoreMs2("", peaks_1$spectra$mz, peaks_1$spectra$intensity,
                    peaks_1$precursor_mz, "", peaks_2$spectra$mz,
                    peaks_2$spectra$intensity, peaks_2$precursor_mz,
                    score_params)
  return(score)

}