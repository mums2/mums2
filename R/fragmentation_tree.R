

#' @export
#' @title compute_molecular_formulas
#' @description
#' Clusters the data together
#' @param mass_data your mass_data object
#' @param parent_ppm ppm
#' @param num_threads data
compute_molecular_formulas <- function(mass_data, parent_ppm = 3, num_threads = detectCores()) {
  size <- length(mass_data$peak_data)
  molecular_formula_list <- vector("list", size)
  pb <- progress_bar$new(
    format = "  Computing Molecular Formulas [:bar] :percent eta: :eta",
    total = size, clear = FALSE, width= 75)
  pb$tick(0)

  for(i in seq_along(mass_data$peak_data)) {
    molecular_formula_list[[i]] <- compute_fragmentation_tree(mass_data$peak_data[[i]],
                                                              mass_data$ms2_matches$mz[[i]],
                                                              parent_ppm, num_threads)
    pb$tick()
  }
  results <- as.character(molecular_formula_list)
  mass_data$predicted_molecular_formulas = results
  return(mass_data)
}


compute_fragmentation_tree <- function(list_of_mz_int, parent_mass, parent_ppm, num_threads) {
  parent_decomp <- decomposeMass(parent_mass, ppm = parent_ppm)
  valid_parent_indexes <- head(which(parent_decomp$valid == "Valid"), 1000)
  invalid_indexes <- head(which(parent_decomp$valid == "Invalid"), 1000)
  if(length(parent_decomp$formula) <= 0) {
    warning("No parent decompositions, returning emptry string.")
    return("")
  }
  if(length(valid_parent_indexes) == 1) {
    return(parent_decomp$formula[valid_parent_indexes])
  }
  if(length(parent_decomp$formula) == 1) {
    warning("No valid parent formulas identified, therefore, return an invalid but possible formula.")
    return(parent_decomp$formula[1])
  }
  # worse case scenario, use the invalid indexes to generate a value
  if(length(valid_parent_indexes) <= 0 && length(invalid_indexes) > 0) {
    warning("No valid parent formulas identified, therefore, replacing with invalid formulas.")
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
