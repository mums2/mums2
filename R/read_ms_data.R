read_mzml_mzxml <- function(file) {
  file_list <- unique(as.list(file))
  all_data <- list(mass_spec_data = data.frame(SpectraIndex = as.integer(),
                                               basePeakMZ = as.numeric(),
                                               retentionTime = as.numeric(),
                                               file = as.character()),
                   peak_data = list())
  for (file in file_list) {
    if (!file.exists(file)) {
      stop(paste0("file: ", file,
                  " does not exist. Please ensure all files exist."))
    }
    extension <- tail(strsplit(file, split = "\\.")[[1]], 1)
    if (!(tolower(extension) %in% c("mzml", "mzxml"))) {
      stop(paste0("Please ensure the input file is a .mzml/mzxml,
                  it is currently a .", extension))
    }

    print(paste0("Reading: ", file, " ..."))
    file_reader <- openMSfile(file)
    peak_length <- length(peaks(file_reader))
    mass_spec_data <- header(file_reader, 1:peak_length)[c("seqNum",
                                                           "msLevel",
                                                           "basePeakMZ",
                                                           "retentionTime")]

    mass_spec_data$basePeakMZ <- as.numeric(mass_spec_data$basePeakMZ)
    mass_spec_data$retentionTime <- as.numeric(mass_spec_data$retentionTime)
    mass_spec_data$file <- rep(file, times = peak_length)
    mass_spec_data$SpectraIndex <- 1:peak_length
    mass_spec_data <- mass_spec_data[mass_spec_data$msLevel == 2, ]

    peak_data <- vector("list", nrow(mass_spec_data))
    for (i in seq_along(peak_data)) {
      peak_data[[i]] <- as.data.frame(peaks(file_reader,
                                            mass_spec_data$seqNum[i]))
    }
    close(file_reader)
    all_data$mass_spec_data <- rbind(all_data$mass_spec_data, mass_spec_data)
    all_data$peak_data[[file]] <- peak_data
  }
  all_data

}

read_mgf <- function(file) {
  file_list <- unique(as.list(file))
  all_data <- list(mass_spec_data = data.frame(SpectraIndex = as.integer(),
                                               basePeakMZ = as.numeric(),
                                               retentionTime = as.numeric(),
                                               file = as.character()),
                   peak_data = list())
  for (file in file_list) {
    if (!file.exists(file)) {
      stop(paste0("file: ", file,
                  " does not exist. Please ensure all files exist."))
    }
    extension <- tail(strsplit(file, split = "\\.")[[1]], 1)
    if (tolower(extension) != "mgf") {
      stop(paste0("Please ensure the input file is a
                  .mgf, it is currently a .", extension))
    }
    print(paste0("Reading: ", file, " ..."))
    result_data_list <- ReadMgf(file)
    row_length <- nrow(result_data_list$ms2_table)
    columns <- colnames(result_data_list$ms2_table)
    rt_type <- which(columns %in% c("RTINMINUTES", "RTINSECONDS"))
    filtered_df <- data.frame(SpectraIndex = 1:row_length,
                              file = rep(file, time = row_length),
                              result_data_list$ms2_table[c("PEPMASS",
                                                           columns[rt_type])])
    colnames(filtered_df) <- c("SpectraIndex", "file",
                               "basePeakMZ", "retentionTime")

    all_data$mass_spec_data <- rbind(all_data$mass_spec_data, filtered_df)
    all_data$peak_data[[file]] <- result_data_list$mzIntensityList
  }
  all_data
}

#' @title Create Reference Database
#' @export
#' @description Reader function msp files
#' @param msp_file the file path of your msp file
#' @examples
#' read_msp(mums2_example("massbank_example_data.msp"))
#'
#' @return a `reference_database` object.
read_msp <- function(msp_file) {
  extension <- tail(strsplit(msp_file, split = "\\.")[[1]], 1)
  if (tolower(extension) != "msp") {
    stop(paste0("Please ensure the input file is a msp,
                it is currently a .", extension))
  }
  print(paste0("Reading: ", msp_file, " ..."))
  reference <- ReadMsp(msp_file)
  class(reference) <- "reference_database"
  return(reference)
}



#' @title Read msp files
#' @export
#' @param hmdb_file the xml hmdb file
#' @param ms2_folder the folder path of your ms2 files
#' @examples
#' read_msp(mums2_example("massbank_example_data.msp"))
#'
#' @return a `reference_database` object.
read_hmdb <- function(hmdb_file, ms2_folder) {
  database <- process_xml(hmdb_file)
  read_and_match_spectra_files(ms2_folder, database)
  annotations <- CreateAnnotationController(database)
  class(annotations) <- "reference_database"
  return(annotations)
}


process_xml <- function(xml_file) {
  print("Reading Metabolites from XML Files...")
  records <- xml_find_all(read_xml(xml_file), "//d1:metabolite")
  print("Processing XML Files...")
  pb <- CreateProgressBarObject()
  database <- CreateHumanMetabolomicsDB()
  progress <- 0
  size <- length(records)
  for(xml in records) {
    tags <- xml_name(xml_children(xml))
    data <- xml_text(xml_children(xml))
    AddHumanMetabolomicNode(database, tags, data)
    IncrementProgressBar(pb, progress / size)
    progress <- progress + 1
  }
  DestroyProgressBar(pb)
  rm(pb)
  return(database)
}


read_and_match_spectra_files <- function(ms2_files, database) {
  ls <- list.files(ms2_files, full.names = T)
  database_names <- sub("_.*", "", list.files(ms2_files, full.names = F))
  AddSpectra(database, ls, database_names)
  ProcessMs2Files(database)
}
