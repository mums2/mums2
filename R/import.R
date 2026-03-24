#' @export
#' @title Import all data
#' @description This function is a wrapper for the mpactr import_data function.
#' It will import your peak table and meta data and create a mpactr_object.
#' @param peak_table The file path to your feature table
#' file.
#' @param meta_data The file path to your meta_data file or `data.frame`.
#' @param format The expected exported type of your peak table, can be
#' one of "Progenesis", "None", "None".
#' @examples
#' data <-
#'    import_all_data(peak_table =
#'                    mums2::mums2_example("botryllus_pt_small.csv"),
#'                    meta_data =
#'                    mums2::mums2_example("meta_data_boryillus.csv"),
#'                    format = "None")
#' @returns a `mpactr` object.
import_all_data <- function(peak_table, meta_data, format) {
  format_to_uft8_remove_commas(import_data(peak_table = peak_table,
                                           meta_data = meta_data,
                                           format = format))
}


#' @export
#' @title Change RT time to minutes or seconds
#' @description This function will change your ms1 peak table rt time to
#' rt time in seconds or minutes. This modification happens in place
#' (or by reference), so the `mpactr_object` will be updated.
#' @param mpactr_object The object created from `import_all_data()`
#' file.
#' @param rt_type how you want to convert your retention time,
#' your options are minutes, or seconds. defaults to seconds.
#' @examples
#' data <-
#'    import_all_data(peak_table =
#'                    mums2::mums2_example("botryllus_pt_small.csv"),
#'                    meta_data =
#'                    mums2::mums2_example("meta_data_boryillus.csv"),
#'                    format = "None")
#' change_rt_to_seconds_or_minute(data, "minutes")
#' @returns a modified `mpactr` object.
change_rt_to_seconds_or_minute <- function(mpactr_object,
                                           rt_type = "seconds") {
  if (!inherits(mpactr_object, "filter_pactr")) {
    stop("Make sure you are using the object created from `import_all_data()`")
  }

  if (rt_type != "seconds" && rt_type != "minutes") {
    stop("rt_type can only be equal to 'seconds', or 'minutes'")
  }
  print(paste0("Changing rt values to ", rt_type))
  peak_table <- get_peak_table(mpactr_object)
  colnames(peak_table)

  if (rt_type == "seconds" && !("RTINSECONDS" %in% colnames(peak_table))) {
    column_index <- c(which(colnames(peak_table) == "rt"),
                      which(colnames(peak_table) == "RTINMINUTES"))
    if (length(column_index) <= 0) {
      stop("There are no colnames that are equal to rt or RTINMINUTES")
    }
    peak_table[[column_index]] <- peak_table[[column_index]] * 60
    colnames(peak_table)[column_index] <- "RTINSECONDS"
    mpactr_object$mpactr_data$set_peak_table(peak_table)
  }
  if (rt_type == "minutes" && !("RTINMINUTES" %in% colnames(peak_table))) {
    column_index <- c(which(colnames(peak_table) == "rt"),
                      which(colnames(peak_table) == "RTINSECONDS"))
    if (length(column_index) <= 0) {
      stop("There are no colnames that are equal to rt or RTINSECONDS")
    }
    peak_table[[column_index]] <- peak_table[[column_index]] / 60
    colnames(peak_table)[column_index] <- "RTINMINUTES"
    mpactr_object$mpactr_data$set_peak_table(peak_table)
  }
  mpactr_object
}



format_to_uft8_remove_commas <- function(mpactr_object) {
  message("If peak table has corrupted compound names they will be converted to
      utf-8 and if there are any commas, they will be converted to periods(.).")
  peak_table <- get_peak_table(mpactr_object)
  # Converts non-utf8 data to utf8 data
  peak_table$Compound <- iconv(peak_table$Compound, from = "latin1", "UTF-8")
  peak_table$Compound <- gsub(",", ".", peak_table$Compound)
  mpactr_object$mpactr_data$set_peak_table(peak_table)
  mpactr_object
}
