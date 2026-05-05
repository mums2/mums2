#' @export
#' @title Compute Molecular formula Other
#' @description
#' de novo algorithm for computing molecular formulas. Using fragmentation trees
#' we are able to generate a resultant molecular formula. To ensure efficient
#' we are using a greedy heurstic to generate the resultant formula. Although
#' this may not always result in the correct prediction, it allows us to
#' efficiently calculate a multitudeof chemical formulas.
#' @param mass_data your mass_data object generated from `ms2_ms1_compare()`
#' @param parent_ppm the ppm you wish to generate the candidate
#'  molecular formulas.
#' @param number_of_threads the amount of threads we the algorithm will use.
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
#' compute_molecular_formulas(matched_data)
#' @references
#' Sebastian Böcker, Florian Rasche, Towards de novo identification of
#' metabolites by analyzing tandem mass spectra, Bioinformatics, Volume 24,
#' Issue 16, August 2008, Pages i49–i55,
#' https://doi.org/10.1093/bioinformatics/btn270
#'
#' @return your mass_data object with an additional `character`
#'  vector of all the predicted formulas.
compute_molecular_formulas <- function(mass_data, parent_ppm = 3,
                                       number_of_threads = detectCores() - 1) {
  if (!inherits(mass_data, "mass_data")) {
    stop(paste0("The mass_data object must be created using the",
                " `ms2_ms1_compare()`"))
  }
  if (!is.numeric(parent_ppm)) {
    stop("parent_ppm must be numeric")
  }

  if (!is.numeric(number_of_threads)) {
    stop("number_of_threads must be numeric")
  }

  if (nrow(mass_data$ms2_matches) <= 0) {
    stop("Your mass_data object has no ms2 matches, cannot continue.")
  }

  if (length(mass_data$peak_data) <= 0) {
    stop("Your mass_data object has no peak data, cannot continue.")
  }
  mzs <- lapply(mass_data$peak_data, function(x) x$mz)
  resultant_formulas <-
    DeNovoMolecularFormulaPrediction(mass_data$ms2_matches$mz, mzs,
                                     parent_ppm, number_of_threads)
  resultant_formulas[which(resultant_formulas == "")] <- NA_character_
  failed_amount <- length(which(is.na(resultant_formulas)))
  message(paste0(abs(length(resultant_formulas) - failed_amount),
                 "/", length(resultant_formulas),
                 " chemical formulas were predicted"))
  mass_data$predicted_molecular_formulas <- resultant_formulas
  mass_data
}
