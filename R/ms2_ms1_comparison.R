#' @title Match your ms1 spectra to a ms2
#' @export
#' @description We are matching your ms1 to your supplied ms2 by looking at the difference between the mz and rt
#' @param ms2_files a list of either *all* mgf files, mzml, or mzxml files. 
#' @param mpactr_object your mpactr object creatd from `import_all_data()`
#' @param mz_tolerance your mass-charge ratio tolerance
#' @param rt_tolerance your retention time tolerance
ms2_ms1_compare <- function(ms2_files, mpactr_object, mz_tolerance, rt_tolerance) {

  ms2_data <- list()
  extension <- tail(strsplit(as.list(ms2_files)[[1]], split = "\\.")[[1]], 1)
  if(tolower(extension) != "mgf") {
    ms2_data <- read_mzml_mzxml(ms2_files)
  }
  else {
    ms2_data <- read_mgf(ms2_files)
  }
  mz2 <- as.numeric(ms2_data$mass_spec_data$basePeakMZ)
  rt2 <- as.numeric(ms2_data$mass_spec_data$retentionTime)
  spectra_index <- ms2_data$mass_spec_data$SpectraIndex
  files <- ms2_data$mass_spec_data$file

  ms1_peak_table <- get_peak_table(mpactr_object)
  mz1 <- ms1_peak_table$mz
  rt1 <- ms1_peak_table$rt
  ms1_compounds <- ms1_peak_table$Compound
  len <- length(ms1_compounds)
  result <- vector("list", len)
  for(i in seq_along(1:len)) {
    mz_error <- abs(mz1[[i]] - mz2) * 1e6 / mz1[[i]] 
    rt_err <- abs(rt1[[i]] - rt2) / rt1[[i]]
    result[[i]] <- which(mz_error <= mz_tolerance & rt_err <= rt_tolerance)
  }
  matched_peaks <- length(which(lapply(result, length) > 0))
  print(paste0(matched_peaks, "/", len, " peaks have an MS2 spectra."))

  resultant_mat <- matrix(0, nrow = matched_peaks, ncol = 5)
  ms2_peaks <- vector("list", matched_peaks)
  colnames(resultant_mat) <- c("mz", "rt", "ms1_compound_id", "spectra_index", "ms2_spectrum_id")
  row_index <- 1
  for(i in seq_along(1:length(result))) {
    if(length(result[[i]]) <= 0)
        next
    index <- result[[i]][which.max(mz2[result[[i]]])]
    mz <- mz2[index]
    rt <- rt2[index]
    resultant_mat[row_index, 1] <- mz
    resultant_mat[row_index, 2] <- rt
    resultant_mat[row_index, 3] <- ms1_compounds[[i]]
    resultant_mat[row_index, 4] <- row_index
    resultant_mat[row_index, 5] <- paste0("mz",mz,"rt",rt)
    ms2_peaks[[row_index]] <- ms2_data$peak_data[[files[index]]][[spectra_index[index]]]
    row_index <- row_index + 1 
  }
  match_df <- as.data.frame(resultant_mat)
  match_df$mz <- as.numeric(match_df$mz)
  match_df$rt <- as.numeric(match_df$rt)
  match_df$spectra_index <- as.numeric(match_df$spectra_index)
  result <- list(ms2_matches = match_df, peak_data = ms2_peaks,
                 ms1_data = ms1_peak_table, samples = get_meta_data(mpactr_object)$Injection)
  class(result) <- "mass_data"
  return(result)
}