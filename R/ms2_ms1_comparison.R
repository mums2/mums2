#' @title Match your ms1 spectra to a ms2
#' @export
#' @description We are matching your ms1 to your supplied ms2 by looking at the difference between the mz and rt
#' @param ms2_data the created object from `read_mzml_mzxml()` or `read_mgf()`
#' @param ms1_peak_table your ms1 peak_table that is created from `import_all_data()`
#' @param mz_tolerance your mass-charge ratio tolerance
#' @param rt_tolerance your retention time tolerance
ms2_ms1_compare <- function(ms2_data, mpactr_object, mz_tolerance, rt_tolerance) {
  mz2 <- as.numeric(ms2_data$mass_spec_data$basePeakMZ)
  rt2 <- as.numeric(ms2_data$mass_spec_data$retentionTime)
  spectra_index <- ms2_data$mass_spec_data$SpectraIndex
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
  colnames(resultant_mat) <- c("mz", "rt", "mz1_compound_id", "spectra_index", "ms2_spectrum_id")
  rowIndex <- 1
  for(i in seq_along(1:length(result))) {
    if(length(result[[i]]) <= 0)
        next
    index <- result[[i]][which.max(mz2[result[[i]]])]
    mz <- mz2[index]
    rt <- rt2[index]
    resultant_mat[rowIndex, 1] <- mz
    resultant_mat[rowIndex, 2] <- rt
    resultant_mat[rowIndex, 3] <- ms1_compounds[[i]]
    resultant_mat[rowIndex, 4] <- spectra_index[index]
    resultant_mat[rowIndex, 5] <- paste0("mz",mz,"rt",rt)
    rowIndex <- rowIndex + 1 
  }
  match_df <- as.data.frame(resultant_mat)
  match_df$mz <- as.numeric(match_df$mz)
  match_df$rt <- as.numeric(match_df$rt)
  match_df$spectra_index <- as.numeric(match_df$spectra_index)
  result <- list(ms2_data = ms2_data, ms2_matches = match_df)
  class(result) <- "mass_data"
  return(result)
}