# Base Filter Function

#' @export
#' @title Filter Peak Table
#' @description This function is a wrapper for all of mpactr's filter functions. 
#' When called with a list of parameters that was generated from one of the following functions, 
#' it will call the subsequent filter: `filter_mispicked_ions_parameters()`, `filter_group_parameters()`,
#' `filter_cv_parameters()`, and `filter_insource_ions_parameters()`. You can also find more information
#' on these functions in`mpactr` documentation.
#' @param mpactr_object the mpactr_object is an object generated from the `import_all_data()` function.
#' This is how we begin our pipeline.
#' @param params the list of arguments generated from calling one of these functions: 
#' `filter_mispicked_ions_parameters()`, `filter_group_parameters()`, `filter_cv_parameters()`,
#' and `filter_insource_ions_parameters()`.
#' @examples
#' data <- import_all_data(peak_table = mums2::mums2_example("full_mix_peak_table.csv"), 
#'                             meta_data = mums2::mums2_example("full_mix_meta_data.csv"), 
#'                              format = "Metaboscape")
#'
#' filtered_data <- data |>
#'    filter_peak_table(filter_mispicked_ions_parameters()) |>
#'    filter_peak_table(filter_cv_parameters(cv_threshold = 0.2)) |>
#'    filter_peak_table(filter_group_parameters(group_threshold = 0.1, "Blanks")) |>
#'    filter_peak_table(filter_insource_ions_parameters())
#' @returns a `mpactr` object that has been filter based on the supplied parameters.
filter_peak_table <- function(mpactr_object, params) {
  return(UseMethod("filter_peak_table", params))
}

# Filter dispatch functions

#' @export
#' @rdname filter_peak_table
filter_peak_table.filter_mispicked_ions <- function(mpactr_object, params) {
  return(filter_mispicked_ions(mpactr_object = mpactr_object,
  ringwin = params$ringwin, isowin = params$isowin,
  trwin = params$trwin, max_iso_shift = params$max_iso_shift,
  merge_peaks = params$merge_peaks, merge_method = params$merge_method,
  copy_object = params$copy_object
  ))
}

#' @export
#' @rdname filter_peak_table
filter_peak_table.filter_group <- function(mpactr_object, params) {
  return(filter_group(mpactr_object = mpactr_object, group_threshold = params$group_threshold,
  group_to_remove = params$group_to_remove, remove_ions = params$remove_ions,
  copy_object = params$copy_object))
}

#' @export
#' @rdname filter_peak_table
filter_peak_table.filter_cv <- function(mpactr_object, params) {
  return(filter_cv(mpactr_object = mpactr_object, cv_threshold = params$cv_threshold,
                  copy_object = params$copy_object))
}

#' @export
#' @rdname filter_peak_table
filter_peak_table.filter_insource_ions <- function(mpactr_object, params) {
  return(filter_insource_ions(mpactr_object = mpactr_object,
                              cluster_threshold = params$cluster_threshold,
                              copy_object = params$copy_object))
}

# Filter Parameter Getters

#' @export
#' @title Filter Mispicked Ions Parameters
#' @description Creates a list of filter mispicked ions arguments
#' for the `filter_peak_table()` function
#' @param ringwin Ringing mass window or detector saturation mass window.
#' Default = 0.5 atomic mass units (AMU).
#' @param isowin Isotopic mass window. Default = 0.01 AMU.
#' @param trwin A `numeric` denoting the retention time threshold for assessing
#' if ions should be merged. Default = 0.005.
#' @param max_iso_shift A `numeric`. Default = 3.
#' @param merge_peaks A `boolean` parameter to determine if peaks found to
#' belong to the same ion should be merged in the feature table.
#' @param merge_method If merge_peaks is TRUE, a method for how similar peaks
#' should be merged. Can be one of "sum".
#' @param copy_object A `boolean` parameter that allows users to return a copied
#' object instead of modifying the object.
#' @examples
#' filter_mispicked_ions_parameters()
#' 
#' @return a `list` object of arguments needed to call the given mpactr function when supplied to the 
#' `filter_peak_table()` wrapper function.
filter_mispicked_ions_parameters <- function(ringwin = 0.5, isowin = 0.01, 
                                             trwin = 0.005, max_iso_shift = 3, 
                                             merge_peaks = TRUE, merge_method = "sum", 
                                             copy_object = FALSE) {
  params <- list(ringwin = ringwin,
                 isowin = isowin, trwin = trwin, max_iso_shift = max_iso_shift,
                 merge_peaks = merge_peaks, merge_method = merge_method,
                 copy_object = copy_object)
  class(params) <- "filter_mispicked_ions"
  return(params)
}

#' @export
#' @title Filter Group Parameters
#' @description Creates a list of filter group arguments
#' for the `filter_peak_table()` function
#' @param group_threshold Relative abundance threshold at which to remove ions.
#' Default = 0.01.
#' @param group_to_remove Biological group name to remove ions from.
#' @param remove_ions A `boolean` parameter. If `TRUE` failing ions will be
#' removed from the peak table. Default = TRUE.
#' @param copy_object A `boolean` parameter that allows users to return a copied
#' object instead of modifying the object.
#' @examples
#' filter_group_parameters(group_to_remove = "blank")
#' @return a `list` object of arguments needed to call the given mpactr function when supplied to the 
#' `filter_peak_table()` wrapper function.
filter_group_parameters <- function(group_threshold = 0.01, group_to_remove,
                                   remove_ions = TRUE, copy_object = FALSE) {
  params <- list(group_threshold = group_threshold,
    group_to_remove = group_to_remove, remove_ions = remove_ions,
    copy_object = copy_object)
  class(params) <- "filter_group"
  return(params)
}

#' @export
#' @title Filter Cv Parameters
#' @description Creates a list of filter cv arguments
#' for the `filter_peak_table()` function
#' @param cv_threshold Coefficient of variation threshold.
#' A lower cv_threshold will result in more stringent filtering and higher
#' reproducibility. Recommended values between 0.2 - 0.5.
#' @param copy_object A `boolean` parameter that allows users to return a copied
#' object instead of modifying the object.
#' @examples
#' filter_cv_parameters(0.2)
#' filter_cv_parameters(0.2)
#' @return a `list` object of arguments needed to call the given mpactr function when supplied to the 
#' `filter_peak_table()` wrapper function.
filter_cv_parameters <- function(cv_threshold = NULL, copy_object = FALSE) {
  params <- list(cv_threshold = cv_threshold, copy_object = copy_object)
  class(params) <- "filter_cv"
  return(params)
}

#' @export
#' @title Filter Insource Ions Parameters
#' @description Creates a list of filter insource ions arguments
#' for the `filter_peak_table()` function
#' @param cluster_threshold Cluster threshold for ion deconvolution.
#' Default = 0.95.
#' @param copy_object A `boolean` parameter that allows users to return
#' a copied object instead of modifying the object.
#' @examples
#' filter_insource_ions_parameters()
#' @return a `list` object of arguments needed to call the given mpactr function when supplied to the 
#' `filter_peak_table()` wrapper function.
filter_insource_ions_parameters <- function(cluster_threshold = 0.95, copy_object = FALSE) {
  params <- list(cluster_threshold = cluster_threshold,
                copy_object = copy_object)
  class(params) <- "filter_insource_ions"
  return(params)
}

