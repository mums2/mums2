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


#' @export
#' @title Change RT time to minutes or seconds
#' @description This function will change your ms1 peak table rt time to rt time in seconds or minutes
#' @param mpactr_object The object created from `import_all_data()`
#' file.
#' @param meta_data The rt time you want to convert to, the options are minutes or seconds.
#' @returns a modified `mpactr` object.
change_rt_to_seconds_or_minutes <- function(mpactr_object, rt_type = "seconds") {
  if(!("filter_pactr" %in% class(mpactr_object))) {
    stop("Make sure you are using the object created from `import_all_data()`")
  }
  print(paste0("Changing rt values to ",rt_type))
  peak_table <- get_peak_table(squid_data)
  if(rt_type == "seconds") {
    peak_table$rt <- peak_table$rt * 60
    mpactr_object$mpactr_data$set_peak_table(peak_table)
  }
  if(rt_type == "minutes") {
    peak_table$rt <- peak_table$rt/60
    mpactr_object$mpactr_data$set_peak_table(peak_table)
  }
  return(mpactr_object)
}
