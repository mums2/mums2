# Base Filter Function

#' @export
#' @title filter peak table
#' @description filter data based on mpactr filters
filter_peak_table <- function(params) {
  return(UseMethod("filter_peak_table"))
}

# Filter dispatch functions
#' @export
#' @title filter mispicked ions wrapper
#' @description filter mispicked ions
filter_peak_table.filter_mispicked_ions <- function(params) {
  return(filter_mispicked_ions(mpactr_object = params$mpactr_object,
  ringwin = params$ringwin, isowin = params$isowin,
  trwin = params$trwin, max_iso_shift = params$max_iso_shift,
  merge_peaks = params$merge_peaks, merge_method = parms$merge_method,
  copy_object = params$copy_object
  ))
}

#' @export
#' @title filter group wrapper
#' @description filter group
filter_peak_table.filter_group <- function(params) {
  return(filter_group(mpactr_object = params$mpactr_object, group_threshold = params$group_threshold,
  group_to_remove = params$group_to_remove, remove_ions = params$remove_ions,
  copy_object = params$copy_object))
}

#' @export
#' @title filter cv wrapper
#' @description filter cv
filter_peak_table.filter_cv <- function(params) {
  return(filter_cv(mpactr_object = params$mpactr_object, cv_threshold = params$cv_threshold,
  cv_param = params$cv_param, copy_object = params$copy_object))
}

#' @export
#' @title filter insource ions wrapper
#' @description filter insource ions
filter_peak_table.filter_insource_ions <- function(params) {
  return(filter_insource_ions(mpactr_object = params$mpactr_object,
                              cluster_threshold = params$cluster_threshold,
                              copy_object = params$copy_object))
}

# Filter Parameter Getters

filter_mispicked_ions_parameters <- function(mpactr_object, ringwin = 0.5, isowin = 0.01, 
                                             trwin = 0.005, max_iso_shift = 3, 
                                             merge_peaks = TRUE, merge_method = "sum", 
                                             copy_object = FALSE) {
  params <- list(mpactr_object = mpactr_object, ringwin = ringwin,
                 isowin = isowin, trwin = trwin, max_iso_shift = max_iso_shift,
                 merge_peaks = merge_peaks, merge_method = merge_method,
                 copy_object = copy_object)
  class(params) <- "filter_mispicked_ions"
  return(params)
}

filter_group_parameters <- function(mpactr_object, group_threshold = 0.01, group_to_remove,
                                   remove_ions = TRUE, copy_object = FALSE) {
  params <- list(mpactr_object = mpactr_object, group_threshold = group_threshold,
    group_to_remove = group_to_remove, remove_ions = remove_ions,
    copy_object = copy_object)
  class(params) <- "filter_group"
  return(params)
}

filter_cv_parameters <- function(mpactr_object, cv_threshold = NULL, cv_param, copy_object = FALSE) {
  params <- list(mpactr_object = mpactr_object, cv_threshold = cv_threshold,
    cv_param = cv_param, copy_object = copy_object)
  
  class(params) <- "filter_cv"
  return(params)
}

filter_insource_ions_parameters <- function(mpactr_object, cluster_threshold = 0.95, copy_object = FALSE) {
  params <- list(mpactr_object = mpactr_object,
                cluster_threshold = cluster_threshold,
                copy_object = copy_object)
  class(params) <- "filter_insource_ions"
  return(params)
}

