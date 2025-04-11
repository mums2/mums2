read_mzml_mzxml <- function(file) {
  file_list <- unique(as.list(file))
  all_data <- list(mass_spec_data = data.frame(SpectraIndex = as.integer(),
                                               basePeakMZ = as.numeric(),
                                               retentionTime = as.numeric(),
                                               file = as.character()), 
                  peak_data = list())
  for(file in file_list) {
    extension <- tail(strsplit(file, split = "\\.")[[1]], 1)
    if(!(tolower(extension) %in% c("mzml", "mzxml"))) {
      stop(paste0("Please ensure the input file is a .mzml/mzxml, it is currently a .", extension))
    }
    
    print(paste0("Reading: ", file, " ..."))
    file_reader <- openMSfile(file)
    peak_length <- length(peaks(file_reader))
    mass_spec_data <- header(file_reader,
                            1:peak_length)[c("seqNum", "msLevel", "basePeakMZ", "retentionTime")]
    
    mass_spec_data$basePeakMZ <- as.numeric(mass_spec_data$basePeakMZ)
    mass_spec_data$retentionTime <- as.numeric(mass_spec_data$retentionTime)
    mass_spec_data$file <- rep(file, times = peak_length)
    mass_spec_data$SpectraIndex <- 1:peak_length
    mass_spec_data <- mass_spec_data[mass_spec_data$msLevel == 2, ]
   
    peak_data <- vector("list", nrow(mass_spec_data))
    for(i in seq_along(peak_data)) {
      peak_data[[i]] <- as.data.frame(peaks(file_reader, mass_spec_data$seqNum[i]))
    }
    close(file_reader)
    all_data$mass_spec_data <- rbind(all_data$mass_spec_data, mass_spec_data)
    all_data$peak_data[[file]] <- peak_data
  }
  return(all_data)

}

read_mgf <- function(file) {
  file_list <- unique(as.list(file))
  all_data <- list(mass_spec_data = data.frame(SpectraIndex = as.integer(),
                                               basePeakMZ = as.numeric(),
                                               retentionTime = as.numeric(),
                                               file = as.character()), 
                  peak_data = list())
  for(file in file_list) {
    extension <- tail(strsplit(file, split = "\\.")[[1]], 1)
    if(tolower(extension) != "mgf") {
      stop(paste0("Please ensure the input file is a .mgf, it is currently a .", extension))
    }
    print(paste0("Reading: ", file, " ..."))
    result_data_list <- ReadMgf(file)
    row_length <- nrow(result_data_list$ms2_table)
    columns <- colnames(result_data_list$ms2_table)
    rt_type <- which(columns %in% c("RTINMINUTES", "RTINSECONDS"))
    filtered_df <- data.frame(SpectraIndex = 1:row_length,
                              file = rep(file, time = row_length),
                              result_data_list$ms2_table[c("PEPMASS", columns[rt_type])])
    colnames(filtered_df) <- c("SpectraIndex", "file", "basePeakMZ", "retentionTime")
   
    all_data$mass_spec_data <- rbind(all_data$mass_spec_data, filtered_df)
    all_data$peak_data[[file]] <- result_data_list$mzIntensityList
  }
  return(all_data)
}

#' @title Read msp files
#' @export
#' @description Reader function msp files
#' @param msp_file the file path of your msp file
read_msp <- function(msp_file) {
  return(ReadMsp(msp_file))
}