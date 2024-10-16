#' Extract MS2 data from multiple classes
#'
#' @details
#' `get_peaks_data` will extract precursor mz and scan (mz, intensity)
#'  information for MS2 data stored in `mass_dataset` and `gnps_ref`
#'  data classes. `mass_dataset` is an S4 object where ms2 metadata
#'  and peaks data are stored in the `ms2_data` slot, while `gnps_ref`
#'  is a list of reference scans in which metadata is stored in
#'  `gnps_ref[[i]]$info` and peaks data are stored in `gnps_ref[[1]]$spec`.
#'
#' Regardless of `data` class type, the position of the desired ms2
#'  scan (`idx`) is used to access a single scan.
#'
#'
#' @param data A data structure containing MS2 data. Currently, this
#'  function supports the classes `mass_dataset` and `gnps_ref`.
#' @param idx A `numeric` value denoting the position of the scan in `data`.
#'
#' @return A `list` of class `peaks_data` with two named vectors:
#'  `precursor_mz` and `spectra`. `precursor_mz` is a `numeric` vector
#'  containing the precursor mz. `spectra` is a `data.frame` with two
#'  columns: `mz` and `intensity`.
#'
#' @export
#'
#' @examples
#' # for class "mass_dataset"
#' \dontrun{
#' dat <- tidyMassDemo
#' peaks <- get_peaks_data(dat, 1)}
#'
#' # for class "gnps_ref"
#' \dontrun{
#'  dat <- psu_msmls
#'  peaks <- get_peaks_data(dat, 1)}
get_peaks_data <- function(data, idx) {
  UseMethod("get_peaks_data", data)
}

#' @method get_peaks_data mass_dataset
#' @export
get_peaks_data.mass_dataset <- function(data, idx) {
  pmz <- data@ms2_data[[1]]@ms2_mz[idx]
  spec <- as.data.frame(data@ms2_data[[1]]@ms2_spectra[[idx]])
  peaks_data <- list("precursor_mz" = pmz,
                     "spectra" = spec)
  class(peaks_data) <- "peaks_data"

  return(peaks_data)
}

#' Create peaks_data from MS2 spectra
#'
#' @details
#' `create_peaks_data` will create an object of class `peaks_data`
#' to store an MS2 scan and associated metadata (*e.g.,* precursor mz).
#'
#'
#' @param spectrum A `data.frame` containing an MS2 spectrum with columns
#' `mz` and `intensity`.
#' @param precursor_mz A `numeric` value denoting the precursor mz for the
#' spectrum.
#'
#' @return A `list` of class `peaks_data` with two named vectors:
#'  `precursor_mz` and `spectra`. `precursor_mz` is a `numeric` vector
#'  containing the precursor mz. `spectra` is a `data.frame` with two
#'  columns: `mz` and `intensity`.
#'
#' @export
#'
#' @examples
#' pmz <- 136.062
#' mz <- c(52.0220, 53.0349, 53.9952, 53.9984, 54.0015, 54.0025, 54.0046,
#'           54.0098, 54.0243, 54.0264, 54.8659, 55.0564, 55.9345, 55.9355,
#'           65.0155, 66.0233, 66.0279, 66.0382, 66.8263, 66.8298, 67.0319,
#'           67.0492, 67.0515, 67.0538, 67.0585, 68.0272, 71.9321, 72.9389,
#'           72.9413)
#' int <- c(8, 276, 15, 6, 7, 753, 21, 15, 9, 25, 29, 11, 6, 22, 44, 33,
#'            5, 3959, 41, 24, 18, 9, 171, 91, 5, 19, 13, 1481, 13)
#' peaks <- create_peaks_data(data.frame("mz" = mz,
#'                                       "intensity" = int),
#'                            pmz)
#'
#' # example from a file
#' peaks <- create_peaks_data(read.csv(example("PSUMSMLS_Adenine.csv")),
#'                            136.0620)
create_peaks_data <- function(spectrum, precursor_mz) {
  pd <- list("precursor_mz" = precursor_mz,
             "spectra" = spectrum)
  class(pd)  <- "peaks_data"

  return(pd)
}

get_ref_precursor <- function(ref_list) {
  pmz_row <- grep("precursorMZ", ref_list$info$key, ignore.case = TRUE)

  return(as.numeric(ref_list$info[pmz_row, 2]))
}
