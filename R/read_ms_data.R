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

    message(paste0("Reading: ", file, " ..."))
    result <- read_ms_data(file)
    all_data$mass_spec_data <- rbind(all_data$mass_spec_data, result$mass_spec_data)
    all_data$peak_data[[file]] <- result$peak_data
  }
  
  all_data
}


read_ms_data <- function(file) {
  mzml_dat <- grabMSdata(file, "MS2")
  mzml_dat$MS2$seq_id <- 1
  counter <- 1
  current_mz <- mzml_dat$MS2$premz[[1]]
  current_rt <- mzml_dat$MS2$rt[[1]]
  for(i in seq(nrow(mzml_dat$MS2))) {
    if(mzml_dat$MS2$premz[[i]] != current_mz ||
      mzml_dat$MS2$rt[[i]] != current_rt) {

      current_mz <- mzml_dat$MS2$premz[[i]]
      current_rt <-  mzml_dat$MS2$rt[[i]]
      counter <- counter + 1
    }
    mzml_dat$MS2$seq_id[[i]] <- counter 
  }

  size <- max(mzml_dat$MS2$seq_id)
  peak_data <- vector("list", size)
  peak_sizes <- table(mzml_dat$MS2$seq_id)
  mass_spec_data <- as.data.frame(matrix(0, nrow = size, ncol = 6))
  colnames(mass_spec_data) <- c("seqNum", "msLevel", "basePeakMZ", "retentionTime", "file", "SpectraIndex")
  idx <- 1
  for(i in seq(size)) {
    # Create Peak Data
    range <- idx:(peak_sizes[[i]] + idx - 1)
    peak_data[[i]]$mz <- mzml_dat$MS2$fragmz[range]
    peak_data[[i]]$intensity <- mzml_dat$MS2$int[range]

    # Create Ms2 Data.frame
    mass_spec_data$seqNum[[i]] <- mzml_dat$MS2$seq_id[[idx]]
    mass_spec_data$msLevel[[i]] <- 2
    mass_spec_data$basePeakMZ[[i]] <- mzml_dat$MS2$premz[[idx]]
    mass_spec_data$retentionTime[[i]] <- mzml_dat$MS2$rt[[idx]]
    mass_spec_data$SpectraIndex[[i]] <- i
    idx <- idx + peak_sizes[[i]]
  }
  mass_spec_data$file <- file
  list(mass_spec_data = mass_spec_data, peak_data = peak_data)
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
    message(paste0("Reading: ", file, " ..."))
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
#' @description Creates a reference database by reading
#' a download msp file. These files can be downloaded from
#' sites like https://systemsomicslab.github.io/compms/msdial/main.html#MSP
#' or https://mona.fiehnlab.ucdavis.edu/downloads
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
  message(paste0("Reading: ", msp_file, " ..."))
  reference <- ReadMsp(msp_file)
  class(reference) <- "reference_database"
  reference
}



#' @title Read HMDB database
#' @export
#' @description
#' This function allows you to create an hmdb database. However
#' you are required to supply an xml hmdb file and a folder path
#' that contains all of the ms2 spectras from the hmdb download
#' page https://www.hmdb.ca/downloads.
#' @param hmdb_file the xml hmdb file.
#' @param ms2_folder the folder path of your ms2 spectra files.
#' @examples
#' read_msp(mums2_example("massbank_example_data.msp" ))
#'
#' @references
#' Wishart DS, Guo A, Oler E, Wang F, Anjum A, Peters H, Dizon R,
#' Sayeeda Z, Tian S, Lee BL, Berjanskii M, Mah R, Yamamoto M,
#' Jovel J, Torres-Calzada C, Hiebert-Giesbrecht M, Lui VW,
#' Varshavi D, Varshavi D, Allen D, Arndt D, Khetarpal N,
#' Sivakumaran A, Harford K, Sanford S, Yee K, Cao X, Budinski Z,
#'  Liigand J, Zhang L, Zheng J, Mandal R, Karu N,
#' Dambrova M, Schiöth HB, Greiner R, Gautam V. HMDB 5.0:
#' the Human Metabolome Database for 2022.
#' Nucleic Acids Res. 2022 Jan 7;50(D1):D622-D631.
#' doi: 10.1093/nar/gkab1062. PMID: 34986597; PMCID: PMC8728138.
#' @return a `reference_database` object.
read_hmdb <- function(hmdb_file, ms2_folder) {
  if (!file.exists(hmdb_file)) {
    stop(paste0("hmdb file",
                " does not exist. Please ensure all files exist."))
  }
  if (!file.exists(ms2_folder)) {
    stop(paste0("ms2_folder",
                " does not exist. Please ensure all files exist."))
  }
  database <- process_xml(hmdb_file)
  read_and_match_spectra_files(ms2_folder, database)
  message("Creating Annotation Database...")
  annotations <- CreateAnnotationController(database)
  class(annotations) <- "reference_database"
  annotations
}

process_xml <- function(xml_file) {
  message("Reading Metabolites from XML Files...")
  records <- xml_find_all(read_xml(xml_file), "//d1:metabolite")
  message("Processing XML Files...")
  pb <- CreateProgressBarObject()
  size <- length(records)
  database <- CreateHumanMetabolomicsDB(size)
  progress <- 0
  for (xml in records) {
    tags <- xml_name(xml_children(xml))
    data <- xml_text(xml_children(xml))
    AddHumanMetabolomicNode(database, tags, data, progress)
    IncrementProgressBar(pb, progress / size)
    progress <- progress + 1
  }
  DestroyProgressBar(pb)
  rm(pb)
  rm(records)
  database
}


read_and_match_spectra_files <- function(ms2_files, database) {
  ls <- list.files(ms2_files, full.names = TRUE)
  database_names <- sub("_.*", "", list.files(ms2_files, full.names = FALSE))
  message("Adding Spectra Files...")
  AddSpectra(database, ls, database_names)
  message("Processing Spectra Files...")
  ProcessMs2Files(database)
}
