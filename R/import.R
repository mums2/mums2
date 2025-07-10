#' @export
#' @title Import all data
#' @description This function is a wrapper for the mpactr import_data function. 
#' It will import your peak table and meta data and create a mpactr_object.
#' @param peak_table The file path to your feature table
#' file.
#' @param meta_data The file path to your meta_data file or `data.frame`.
#' @param format The expected exported type of your peak table, can be
#' one of "Progenesis", "Metaboscape", "None".
#' @examples
#' data <- import_all_data(peak_table = mums2::mums2_example("full_mix_peak_table.csv"), 
#'                             meta_data = mums2::mums2_example("full_mix_meta_data.csv"), 
#'                              format = "Metaboscape")
#' @returns a `mpactr` object.
import_all_data <- function(peak_table, meta_data, format) {
  return(import_data(peak_table = peak_table,
                     meta_data = meta_data, format = format))
}


#' @export
#' @title Change RT time to minutes or seconds
#' @description This function will change your ms1 peak table rt time to 
#' rt time in seconds or minutes. This modification happens in place (or by reference),
#' so the `mpactr_object` will be updated.
#' @param mpactr_object The object created from `import_all_data()`
#' file.
#' @param rt_type how you want to convert your retention time, your options are minutes, or seconds.
#' defaults to seconds.
#' @examples
#' data <- import_all_data(peak_table = mums2::mums2_example("full_mix_peak_table.csv"), 
#'                             meta_data = mums2::mums2_example("full_mix_meta_data.csv"), 
#'                              format = "Metaboscape")
#' change_rt_to_seconds_or_minutes(data, "minutes")
#' @returns a modified `mpactr` object.
change_rt_to_seconds_or_minutes <- function(mpactr_object, rt_type = "seconds") {
  if(!("filter_pactr" %in% class(mpactr_object))) {
    stop("Make sure you are using the object created from `import_all_data()`")
  }
  print(paste0("Changing rt values to ",rt_type))
  peak_table <- get_peak_table(mpactr_object)
  colnames(peak_table)

  if(rt_type == "seconds" && !("RTINSECONDS" %in% colnames(peak_table))) {
   
    column_index <- c(which(colnames(peak_table) == "rt"), which(colnames(peak_table) == "RTINMINUTES"))
    if(length(column_index) <= 0) {
      stop("There are no colnames that are equal to rt or RTINMINUTES")
    }
    peak_table[[column_index]] <- peak_table[[column_index]] * 60
    colnames(peak_table)[column_index] = "RTINSECONDS"
    mpactr_object$mpactr_data$set_peak_table(peak_table)
  }
  if(rt_type == "minutes" && !("RTINMINUTES" %in% colnames(peak_table))) {
    column_index <- c(which(colnames(peak_table) == "rt"), which(colnames(peak_table) == "RTINSECONDS"))
    if(length(column_index) <= 0) {
      stop("There are no colnames that are equal to rt or RTINSECONDS")
    }
    peak_table[[column_index]] <- peak_table[[column_index]] / 60
    colnames(peak_table)[column_index] = "RTINMINUTES"
    mpactr_object$mpactr_data$set_peak_table(peak_table)
  }
  return(mpactr_object)
}
