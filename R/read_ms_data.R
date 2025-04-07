read_ms_data <- function(file) {
  extension <- tail(strsplit(file, split = "\\.")[[1]], 1)
  if(extension == "mgf") {
    return(read_mgf(file))
  }
  file_reader <- mzR::openMSfile(file)
  mass_spec_data <- header(file_reader, 
                           1:length(peaks(file_reader)))[c("seqNum", "msLevel", "basePeakMZ", "retentionTime")]
  mass_spec_data <- mass_spec_data[mass_spec_data$msLevel == 2, ]
  peak_data <- vector("list", nrow(mass_spec_data))
  for(i in seq_along(peak_data)) {
    peak_data[[i]] <- as.data.frame(peaks(file_reader, mass_spec_data$seqNum[i]))
  }
  mzR::close(file_reader)
  return(list(mass_spec_data = data.frame(SpectraIndex = 1:length(nrow(mass_spec_data)), 
                                          mass_spec_data[c("basePeakMZ", "retentionTime")]),
              peak_data = peak_data))
}

read_mgf <- function(file) {
  result_data_list <- Read(file)
  filtered_df <- data.frame(SpectraIndex = 1:nrow(result_data_list$ms2_table),
                            result_data_list$ms2_table[c("PEPMASS", "RTINMINUTES")])
  colnames(filtered_df) <- c("SpectraIndex", "basePeakMZ", "retentionTime")
  return(list(mass_spec_data = filtered_df, peak_data = result_data_list$mzIntensityList)) 
}
