#' GNPS-like similarity between two MS/MS spectra
#'
#' @description
#' `modified_cosine_params()` generates a parameter list to perform GNPS-like
#' cosine similarity score calculation between two MS2 spectra.
#'
#' @details
#' `modified_cosine_params()` will initiate cosine scoring based on the Python
#' code by Wang et al. (2016), which is currently used for cosine scoring
#' in GNPS, to calculate similarity between two MS2 spectra. This scoring
#' method will compare peaks data, apply a square root normalization
#' to peak intensities, align peaks both with and without correction
#' for mass shifts, and calculate similarity.
#'
#'
#' @param frag_tolerance The mz fragment tolerance threshold for aligning
#' fragment peaks from two ms2 spectra. GNPS default = 0.5.
#'
#' @examples
#' modified_cosine_params(0.5)
#'
#' @return A parameters list for similarity scoring method "gnps"
#' @references
#' Mingxun Wang, Jeremy J. Carver, Vanessa V. Phelan, Laura M. Sanchez,
#' Neha Garg, Yao Peng, Don Duy Nguyen et al. "Sharing and community curation
#' of mass spectrometry data with Global Natural Products Social Molecular
#' Networking." Nature biotechnology 34, no. 8 (2016): 828. PMID: 27504778
#'
#' @export
modified_cosine_params <- function(frag_tolerance) {
  parameters <- list("tolerance" = frag_tolerance,
                     "method" = "gnps")
  class(parameters) <- "parameters"
  parameters
}

#' Entropy similarity between two MS/MS spectra
#'
#' @description
#' Calculate spectral entropy similarity between two MS2 spectra
#'
#'
#' @details
#' `spec_entropy_params()` will initiate spectral entropy similarity scoring via
#' the `msentropy` package (Li et al. 2021). For more information about
#' parameters, see [msentropy::msentropy_similarity()].
#'
#' @param ms2_tolerance_in_da MS2 peak tolerance in Da, set to -1 to disable.
#' Defaults to `0.02`.
#' @param ms2_tolerance_in_ppm MS2 peak tolerance in ppm, set to -1 to disable.
#' Defaults to `-1`.
#' @param clean_spectra Either `TRUE` or `FALSE` to clean the spectra prior to
#' calculating similarity. See `msentropy::clean_spectrum` for more information.
#' Defaults to `TRUE`.
#' @param min_mz `numeric`, minimum mz to keep, set to -1 to disable. Defaults
#' to `0`.
#' @param max_mz `numeric`, maximum mz to keep, set to -1 to disable. Defaults
#' to `1000`.
#' @param noise_threshold Background intensity threshold, all peaks with
#' intensity < noise_threshold * max_intensity are removed. Set to -1 to
#' disable. Defaults to `0.01`.
#' @param max_peak_num `numeric`, maximum number of peaks to keep for score
#' calculation. Set to -1 to disable. Defaults to `100`.
#' @param weighted `logical` whether weighted or unweighted entropy similarity
#' will be calculated. Defaults to `TRUE`.
#' @examples
#' spec_entropy_params()
#' @return A parameters `list` for similarity scoring method "spectral_entropy"
#'
#' @export
#'
#' @references
#' Li, Y., Kind, T., Folz, J. et al. Spectral entropy outperforms MS/MS dot
#' product similarity for small-molecule compound identification, Nat Methods
#' 18, 1524–1531 (2021). https://doi.org/10.1038/s41592-021-01331-z
spec_entropy_params <- function(ms2_tolerance_in_da = 0.02,
                                ms2_tolerance_in_ppm = -1,
                                clean_spectra = TRUE, min_mz = 0, max_mz = 1000,
                                noise_threshold = 0.01, max_peak_num = 100,
                                weighted = TRUE) {
  parameters <- list("ms2_tolerance_in_da" = ms2_tolerance_in_da,
                     "ms2_tolerance_in_ppm" = ms2_tolerance_in_ppm,
                     "clean_spectra" = clean_spectra,
                     "min_mz" = min_mz,
                     "max_mz" = max_mz,
                     "noise_threshold" = noise_threshold,
                     "max_peak_num" = max_peak_num,
                     "weighted" = weighted,
                     "method" = "entropy")
  class(parameters) <- "parameters"
  parameters
}
