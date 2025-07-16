

#' @export
#' @title Compute Molecular formula
#' @description
#' De-novo algorithm for computing molecular formulas. Using fragmentation trees we are able to generate
#' a resultant molecular formula. To ensure efficent we are using a greedy heurstic to generate the resultant formula.
#' Although this may not always result in the correct prediction, it allows us to efficently calculate a multitude
#' of chemical formulas.
#' @param mass_data your mass_data object generated from `ms2_ms1_compare()`
#' @param parent_ppm the ppm you wish to generate the candidate molecular formulas.
#' @param num_threads the amount of threads we the algorithm will use. 
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
#'
#'
#' matched_data <- ms2_ms1_compare(mums2_example("full_mix_ms2.mgf"),
#'  filtered_data, 2, 6)
#' compute_molecular_formulas(matched_data)
#' @return your mass_data object with an additional `character` vector of all the predicted formulas.
compute_molecular_formulas <- function(mass_data, parent_ppm = 3, num_threads = detectCores()) {
  size <- length(mass_data$peak_data)
  molecular_formula_list <- vector("list", size)
  pb <- CreateProgressBarObject()

  for(i in seq_along(mass_data$peak_data)) {
    molecular_formula_list[[i]] <- compute_fragmentation_tree(mass_data$peak_data[[i]],
                                                              mass_data$ms2_matches$mz[[i]],
                                                              parent_ppm, num_threads)
    IncrementProgressBar(pb, i/size)
  }
  DestroyProgressBar(pb)
  rm(pb)
  results <- as.character(molecular_formula_list)
  mass_data$predicted_molecular_formulas = results
  failed_amount <- length(which(is.na(results)))
  message(paste0(abs(length(results) - failed_amount), "/", length(results),
                 " chemical formulas were predicted"))
  return(mass_data)
}


compute_fragmentation_tree <- function(list_of_mz_int, parent_mass, parent_ppm, num_threads) {
  parent_decomp <- decomposeMass(parent_mass, ppm = parent_ppm)
  valid_parent_indexes <- head(which(parent_decomp$valid == "Valid"), 1000)
  invalid_indexes <- head(which(parent_decomp$valid == "Invalid"), 1000)
  if(length(parent_decomp$formula) <= 0) {
    return(NA_character_)
  }
  if(length(valid_parent_indexes) == 1) {
    return(parent_decomp$formula[valid_parent_indexes])
  }
  if(length(parent_decomp$formula) == 1) {
    return(parent_decomp$formula[1])
  }
  # worse case scenario, use the invalid indexes to generate a value
  if(length(valid_parent_indexes) <= 0 && length(invalid_indexes) > 0) {
    valid_parent_indexes <- invalid_indexes
  }
  decompList <- vector("list", length(list_of_mz_int$mz))
  for(i in seq_along(list_of_mz_int$mz)){
    if(parent_mass < list_of_mz_int$mz[[i]]){
      next
    }
    decompList[[i]] <- decomposeIsotopes(list_of_mz_int$mz[[i]], list_of_mz_int$intensity[[i]])
  }
  full_data <- c()
  color_count <- 1
  for(i in seq_along(decompList)){
    if(is.null(decompList[[i]])) {
      next
    }
    valid_indexes <- head(which(decompList[[i]]$valid == "Valid"), 1000)
    scores <- decompList[[i]]$score[valid_indexes]
    full_data$score <- append(full_data$score, scores)
    full_data$formula <- append(full_data$formula, decompList[[i]]$formula[valid_indexes])
    full_data$mass <- append(full_data$mass, decompList[[i]]$exactmass[valid_indexes])
    full_data$color <- c(full_data$color, rep(color_count, length(scores)))
    color_count <- color_count + 1
  }
  if(length(full_data$score) <= 0) {
    return(parent_decomp$formula[[1]])
  }
  scores <- parent_decomp$score[valid_parent_indexes]
  full_data$score <- append(scores, full_data$score)
  full_data$formula <- append(parent_decomp$formula[valid_parent_indexes], full_data$formula)
  full_data$mass <- append(parent_decomp$exactmass[valid_parent_indexes], full_data$mass)
  full_data$color <- append(rep(0, length(scores)), full_data$color)
  res <- ComputeFragmentationTree(full_data, parent_mass, num_threads)
  return(res)
}
