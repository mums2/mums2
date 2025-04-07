#' @title Read mzml and mzXML files
#' @export
#' @description Reader function for mzml and mzXML files
#' @param file the object path of your mzml/mzXML file
read_mzml_mzxml <- function(file) {
  extension <- tail(strsplit(file, split = "\\.")[[1]], 1)
  if(!(tolower(extension) %in% c("mzml", "mzxml"))) {
    stop(paste0("Please ensure the input file is a .mzml/mzxml, it is currently a .", extension))
  }
  file_reader <- openMSfile(file)
  mass_spec_data <- header(file_reader, 
                           1:length(peaks(file_reader)))[c("seqNum", "msLevel", "basePeakMZ", "retentionTime")]
  mass_spec_data$basePeakMZ <- as.numeric(mass_spec_data$basePeakMZ)
  mass_spec_data$retentionTime <- as.numeric(mass_spec_data$retentionTime)
  
  mass_spec_data <- mass_spec_data[mass_spec_data$msLevel == 2, ]
  peak_data <- vector("list", nrow(mass_spec_data))
  for(i in seq_along(peak_data)) {
    peak_data[[i]] <- as.data.frame(peaks(file_reader, mass_spec_data$seqNum[i]))
  }
  close(file_reader)
  return(list(mass_spec_data = data.frame(SpectraIndex = 1:nrow(mass_spec_data), 
                                          mass_spec_data[c("basePeakMZ", "retentionTime")]),
              peak_data = peak_data))
}

#' @title Read mgf files
#' @export
#' @description Reader function mgf files
#' @param file the object path of your mgf file
read_mgf <- function(file) {
  extension <- tail(strsplit(file, split = "\\.")[[1]], 1)
  if(tolower(extension) != "mgf") {
    stop(paste0("Please ensure the input file is a .mgf, it is currently a .", extension))
  }
  result_data_list <- Read(file)
  filtered_df <- data.frame(SpectraIndex = 1:nrow(result_data_list$ms2_table),
                            result_data_list$ms2_table[c("PEPMASS", "RTINMINUTES")])
  colnames(filtered_df) <- c("SpectraIndex", "basePeakMZ", "retentionTime")
  return(list(mass_spec_data = filtered_df, peak_data = result_data_list$mzIntensityList)) 
}
