#' @export
#' @title Import all data
#' @description This function is a wrapper for the mpactr import_data function. It will import your peak table and meta data and create a mpactr_object.
#' @param peak_table The file path to your feature table
#' file.
#' @param meta_data The file path to your meta_data file or `data.frame`.
#' @param format The expected exported type of your peak table, can be
#' one of "Progenesis", "Metaboscape", "None".
#' @returns a `mpactr` object.
import_all_data <- function(peak_table, meta_data, format) {
  return(import_data(peak_table = peak_table,
                     meta_data = meta_data, format = format))
}
