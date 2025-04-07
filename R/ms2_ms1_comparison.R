
ms2_ms1_compare <- function(ms2_data, ms1_data, rt_tolerance, mz_tolerance) {
  mz_tolerance <- 100
  rt_tolerance <- 150
  ms2_data <- dat$mass_spec_data
  ms1_data <- ms1_peaks
  mz2 <- as.numeric(ms2_data$basePeakMZ)
  rt2 <- as.numeric(ms2_data$retentionTime)

  mz1 <- ms1_data$mz
  rt1 <- ms1_data$rt
  ms1_compounds <- ms1_data$Compound

  len <- nrow(mass_data_set@variable_info)
  result <- vector("list", len)
  for(i in seq_along(1:len)) {
    mz_error <- abs(mz1[[i]] - mz2) * 1e6 / mz1[[i]] 
    rt_err <- abs(rt1[[i]] - rt2) / rt1[[i]]
    result[[i]] <- which(mz_error <= mz_tolerance & rt_err <= rt_tolerance)
  }
  matched_peaks <- length(which(lapply(result, length) > 0))
  print(paste0(matched_peaks, "/", len, " peaks have an MS2 spectra."))

  resultant_mat <- matrix(0, nrow = matched_peaks, ncol = 3)
  colnames(resultant_mat) <- c("mz", "rt", "mz1_compound_id")
  rowIndex <- 1
  for(i in seq_along(1:length(result))) {
    if(length(result[[i]]) <= 0)
        next
    index <- which.max(mz2[result[[i]]])
    resultant_mat[rowIndex, 1] <- mz2[result[[i]][[index]]]
    resultant_mat[rowIndex, 2] <- rt2[result[[i]][[index]]]
    resultant_mat[rowIndex, 3] <- ms1_compounds[[i]]
    rowIndex <- rowIndex + 1
  }
  return(resultant_mat)
}

microbenchmark::microbenchmark(abs(mz1[[2]] - mz2), test(mz2, mz1[[2]]))
test(mz2, mz1[[1]])

#' @export
test <- function(mz2, mz1)
{
  Diffs(mz2, mz1)
}

vals <- length(which(lapply(result, length) > 0))
all(resultant_mat[,3] == mass_data_set@ms2_data$`12152023_Coculture_with_new_JC1.gnps.mgf`@variable_id)
