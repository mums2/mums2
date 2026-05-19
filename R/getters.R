#' @export
#' @title Get MS2 Matches
#' @description
#' A getter that returns the generated ms2 data in your `mass_data` object.
#' @param mass_data The object generated from `ms2_ms1_compare()`.
#' @examples
#' data <-
#'    import_all_data(peak_table =
#'                    mums2::mums2_example("botryllus_pt_small.csv"),
#'                    meta_data =
#'                    mums2::mums2_example("meta_data_boryillus.csv"),
#'                    format = "None")
#'
#' filtered_data <- data |>
#'    filter_peak_table(filter_mispicked_ions_params()) |>
#'    filter_peak_table(filter_cv_params(cv_threshold = 0.2)) |>
#'    filter_peak_table(filter_group_params(group_threshold = 0.1,
#'                                              "Blanks")) |>
#'    filter_peak_table(filter_insource_ions_params())
#'
#' change_rt_to_seconds_or_minute(filtered_data, "minutes")
#' 
#' get_ms2_matches(matched_data)
#'
#' @return a `data.frame` object containing ms2 data.
get_ms2_matches <- function(mass_data) {
  if(!inherits(mass_data, "mass_data")) {
    stop("mass_data must be generated from the `ms2_ms1_compare()` function.")
  }
  mass_data$ms2_matches
}

#' @export
#' @title Get MS1 Matches
#' @description
#' A getter that returns the generated ms2 data in your `mass_data` object.
#' @param mass_data The object generated from `ms2_ms1_compare()`.
#' @examples
#' data <-
#'    import_all_data(peak_table =
#'                    mums2::mums2_example("botryllus_pt_small.csv"),
#'                    meta_data =
#'                    mums2::mums2_example("meta_data_boryillus.csv"),
#'                    format = "None")
#'
#' filtered_data <- data |>
#'    filter_peak_table(filter_mispicked_ions_params()) |>
#'    filter_peak_table(filter_cv_params(cv_threshold = 0.2)) |>
#'    filter_peak_table(filter_group_params(group_threshold = 0.1,
#'                                              "Blanks")) |>
#'    filter_peak_table(filter_insource_ions_params())
#'
#' change_rt_to_seconds_or_minute(filtered_data, "minutes")
#'
#' matched_data <- ms2_ms1_compare(mums2_example("botryllus_v2.gnps.mgf"),
#'  filtered_data, 10, 6)
#' 
#' get_ms1_data(matched_data)
#'
#' @return a `data.frame` object containing ms1 data.
get_ms1_data <- function(mass_data) {
  if(!inherits(mass_data, "mass_data")) {
    stop("mass_data must be generated from the `ms2_ms1_compare()` function.")
  }
  mass_data$ms1_data
}

#' @export
#' @title Get Peaks Data
#' @description
#' A getter that will return the peak data of all of the matched specturms.
#' The peak data is the list of mz/intensities found in the ms2 file.
#' @param mass_data The object generated from `ms2_ms1_compare()`.
#' @examples
#' data <-
#'    import_all_data(peak_table =
#'                    mums2::mums2_example("botryllus_pt_small.csv"),
#'                    meta_data =
#'                    mums2::mums2_example("meta_data_boryillus.csv"),
#'                    format = "None")
#'
#' filtered_data <- data |>
#'    filter_peak_table(filter_mispicked_ions_params()) |>
#'    filter_peak_table(filter_cv_params(cv_threshold = 0.2)) |>
#'    filter_peak_table(filter_group_params(group_threshold = 0.1,
#'                                              "Blanks")) |>
#'    filter_peak_table(filter_insource_ions_params())
#'
#' change_rt_to_seconds_or_minute(filtered_data, "minutes")
#'
#' matched_data <- ms2_ms1_compare(mums2_example("botryllus_v2.gnps.mgf"),
#'  filtered_data, 10, 6)
#' 
#' get_ms2_peaks_data(matched_data)
#'
#' @return a `list` object containing peak data.
get_ms2_peaks_data <- function(mass_data) {
  if(!inherits(mass_data, "mass_data")) {
    stop("mass_data must be generated from the `ms2_ms1_compare()` function.")
  }
  mass_data$peak_data
}

#' @export
#' @title Get Samples
#' @description
#' Returns a list of your samples found in the metadata file.
#' @param mass_data The object generated from `ms2_ms1_compare()`.
#' @examples
#' data <-
#'    import_all_data(peak_table =
#'                    mums2::mums2_example("botryllus_pt_small.csv"),
#'                    meta_data =
#'                    mums2::mums2_example("meta_data_boryillus.csv"),
#'                    format = "None")
#'
#' filtered_data <- data |>
#'    filter_peak_table(filter_mispicked_ions_params()) |>
#'    filter_peak_table(filter_cv_params(cv_threshold = 0.2)) |>
#'    filter_peak_table(filter_group_params(group_threshold = 0.1,
#'                                              "Blanks")) |>
#'    filter_peak_table(filter_insource_ions_params())
#'
#' change_rt_to_seconds_or_minute(filtered_data, "minutes")
#'
#' matched_data <- ms2_ms1_compare(mums2_example("botryllus_v2.gnps.mgf"),
#'  filtered_data, 10, 6)
#' 
#' get_samples(matched_data)
#'
#' @return a `character` vector contain all of your samples.
get_samples <- function(mass_data) {
  if(!inherits(mass_data, "mass_data")) {
    stop("mass_data must be generated from the `ms2_ms1_compare()` function.")
  }
  mass_data$samples
}

#' @export
#' @title Get Molecular Formula Predictions
#' @description
#' Returns all of the molecular formula predictions.
#' @param mass_data The object generated from `ms2_ms1_compare()`.
#' @examples
#' data <-
#'    import_all_data(peak_table =
#'                    mums2::mums2_example("botryllus_pt_small.csv"),
#'                    meta_data =
#'                    mums2::mums2_example("meta_data_boryillus.csv"),
#'                    format = "None")
#'
#' filtered_data <- data |>
#'    filter_peak_table(filter_mispicked_ions_params()) |>
#'    filter_peak_table(filter_cv_params(cv_threshold = 0.2)) |>
#'    filter_peak_table(filter_group_params(group_threshold = 0.1,
#'                                              "Blanks")) |>
#'    filter_peak_table(filter_insource_ions_params())
#'
#'
#' matched_data <- ms2_ms1_compare(mums2_example("botryllus_v2.gnps.mgf"),
#'  filtered_data, 0.1, 6)
#' matched_data <- compute_molecular_formulas(matched_data)
#' 
#' get_molecular_formula_predictions(matched_data)
#'
#' @return a `character` vector contain all of your predicted molecular formulas.
get_molecular_formula_predictions <- function(mass_data) {
  if(!inherits(mass_data, "mass_data")) {
    stop("mass_data must be generated from the `ms2_ms1_compare()` function.")
  }
  if(is.null(mass_data$predicted_molecular_formulas)) {
    stop(paste0("No data found, make sure you run the", 
    " `compute_molecular_formulas()` function before you call",
    " this function."))
  }
  mass_data$predicted_molecular_formulas
}
