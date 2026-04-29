

#' @export
#' @title Compute Molecular formula
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
                                       number_of_threads = detectCores()) {
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

  size <- length(mass_data$peak_data)
  molecular_formula_list <- vector("list", size)
  pb <- CreateProgressBarObject()

  for (i in seq_along(mass_data$peak_data)) {
    molecular_formula_list[[i]] <-
      compute_fragmentation_tree(mass_data$peak_data[[i]],
                                 mass_data$ms2_matches$mz[[i]],
                                 parent_ppm, number_of_threads)
    IncrementProgressBar(pb, i / size)
  }
  DestroyProgressBar(pb)
  rm(pb)
  results <- as.character(molecular_formula_list)
  mass_data$predicted_molecular_formulas <- results
  failed_amount <- length(which(is.na(results)))
  message(paste0(abs(length(results) - failed_amount), "/", length(results),
                 " chemical formulas were predicted"))
  return(mass_data)
}


compute_fragmentation_tree <- function(list_of_mz_int, parent_mass,
                                       parent_ppm, num_threads) {
  parent_decomp <- decomposeMass(parent_mass, ppm = parent_ppm)
  valid_parent_indexes <- head(which(parent_decomp$valid == "Valid"), 1000)
  invalid_indexes <- head(which(parent_decomp$valid == "Invalid"), 1000)
  if (length(parent_decomp$formula) <= 0) {
    return(NA_character_)
  }
  if (length(valid_parent_indexes) == 1) {
    return(parent_decomp$formula[valid_parent_indexes])
  }
  if (length(parent_decomp$formula) == 1) {
    return(parent_decomp$formula[1])
  }
  # worse case scenario, use the invalid indexes to generate a value
  if (length(valid_parent_indexes) <= 0 && length(invalid_indexes) > 0) {
    valid_parent_indexes <- invalid_indexes
  }
  decomp_list <- vector("list", length(list_of_mz_int$mz))
  for (i in seq_along(list_of_mz_int$mz)) {
    if (parent_mass < list_of_mz_int$mz[[i]]) {
      next
    }
    decomp_list[[i]] <- decomposeIsotopes(list_of_mz_int$mz[[i]],
                                          list_of_mz_int$intensity[[i]])
  }
  full_data <- c()
  color_count <- 1
  for (i in seq_along(decomp_list)){
    if (is.null(decomp_list[[i]])) {
      next
    }
    valid_indexes <- head(which(decomp_list[[i]]$valid == "Valid"), 1000)
    scores <- decomp_list[[i]]$score[valid_indexes]
    full_data$score <- append(full_data$score, scores)
    full_data$formula <- append(full_data$formula,
                                decomp_list[[i]]$formula[valid_indexes])
    full_data$mass <- append(full_data$mass,
                             decomp_list[[i]]$exactmass[valid_indexes])
    full_data$color <- c(full_data$color, rep(color_count, length(scores)))
    color_count <- color_count + 1
  }
  if (length(full_data$score) <= 0) {
    return(parent_decomp$formula[[1]])
  }
  scores <- parent_decomp$score[valid_parent_indexes]
  full_data$score <- append(scores, full_data$score)
  full_data$formula <- append(parent_decomp$formula[valid_parent_indexes],
                              full_data$formula)
  full_data$mass <- append(parent_decomp$exactmass[valid_parent_indexes],
                           full_data$mass)
  full_data$color <- append(rep(0, length(scores)), full_data$color)
  res <- ComputeFragmentationTree(full_data, parent_mass, num_threads)
  res
}






#' @export
#' @title Compute Molecular formula2
#' @description
#' de novo algorithm for computing molecular formulas. Using fragmentation trees
#' we are able to generate a resultant molecular formula. To ensure efficient
#' we are using a greedy heurstic to generate the resultant formula. Although
#' this may not always result in the correct prediction, it allows us to
#' efficiently calculate a multitudeof chemical formulas.
#' @param matched_data your mass_data object generated from `ms2_ms1_compare()`
#' @param parent_ppm the ppm you wish to generate the candidate
#'  molecular formulas.
#' @param number_of_threads the amount of threads we the algorithm will use.
compute_molecular_formulas2 <- function(matched_data, parent_ppm = 3,
                                       number_of_threads = detectCores() - 1) {

  size <- nrow(matched_data$ms2_matches)
  pboptions(use_lb = TRUE, nout = 10)
  cl <- parallel::makePSOCKcluster(getOption("cl.cores", number_of_threads))

  clusterExport(cl, c("matched_data"), envir = environment())
  print("Generating Potential formulas...")
  potential_formulas <- pblapply(seq(size), function(i) {
    mums2:::create_all_possible_formulas(matched_data$peak_data[[i]],
                                  matched_data$ms2_matches$mz[[i]], parent_ppm, i)
  } , cl = cl)
  stopCluster(cl)
  pboptions(use_lb = FALSE, nout = 100)
  molecular_formula_list <- vector("list", size)
  pb <- CreateProgressBarObject()
  for (i in seq_along(potential_formulas)) {
    formulas <- potential_formulas[[i]]
    if(inherits(formulas$full_data, "character") || is.null(formulas$full_data)) {
      molecular_formula_list[formulas$index] <-  formulas$full_data
      next
    }
    molecular_formula_list[[formulas$index]] <-
      ComputeFragmentationTree(formulas$full_data, formulas$parent_mass,
                               number_of_threads)
    IncrementProgressBar(pb, i / size)
  }
  DestroyProgressBar(pb)
  rm(pb)
  results <- as.character(molecular_formula_list)
  matched_data$predicted_molecular_formulas <- results
  failed_amount <- length(which(is.na(results)))
  message(paste0(abs(length(results) - failed_amount), "/", length(results),
                 " chemical formulas were predicted"))
  matched_data
}


create_all_possible_formulas <- function(list_of_mz_int, parent_mass,
                                       parent_ppm, index) {
  parent_decomp <- decomposeMass(parent_mass, ppm = parent_ppm)
  valid_parent_indexes <- head(which(parent_decomp$valid == "Valid"), 1000)
  invalid_indexes <- head(which(parent_decomp$valid == "Invalid"), 1000)
  if (length(parent_decomp$formula) <= 0) {
    return(list(full_data = NA_character_, index = index))
  }
  if (length(valid_parent_indexes) == 1) {
    return(list(full_data = parent_decomp$formula[valid_parent_indexes],
                index = index))
  }
  if (length(parent_decomp$formula) == 1) {
    return(list(full_data = parent_decomp$formula[1], index = index))
  }
  # worse case scenario, use the invalid indexes to generate a value
  if (length(valid_parent_indexes) <= 0 && length(invalid_indexes) > 0) {
    valid_parent_indexes <- invalid_indexes
  }
  decomp_list <- vector("list", length(list_of_mz_int$mz))
  for (i in seq_along(list_of_mz_int$mz)) {
    if (parent_mass < list_of_mz_int$mz[[i]]) {
      next
    }
    decomp_list[[i]] <- decomposeIsotopes(list_of_mz_int$mz[[i]],
                                          list_of_mz_int$intensity[[i]])
  }
  full_data <- c()
  color_count <- 1
  for (i in seq_along(decomp_list)){
    if (is.null(decomp_list[[i]])) {
      next
    }
    valid_indexes <- head(which(decomp_list[[i]]$valid == "Valid"), 1000)
    scores <- decomp_list[[i]]$score[valid_indexes]
    full_data$score <- append(full_data$score, scores)
    full_data$formula <- append(full_data$formula,
                                decomp_list[[i]]$formula[valid_indexes])
    full_data$mass <- append(full_data$mass,
                             decomp_list[[i]]$exactmass[valid_indexes])
    full_data$color <- c(full_data$color, rep(color_count, length(scores)))
    color_count <- color_count + 1
  }
  if (length(full_data$score) <= 0) {
    return(list(full_data = parent_decomp$formula[[1]],
                index = index))
  }
  scores <- parent_decomp$score[valid_parent_indexes]
  full_data$score <- append(scores, full_data$score)
  full_data$formula <- append(parent_decomp$formula[valid_parent_indexes],
                              full_data$formula)
  full_data$mass <- append(parent_decomp$exactmass[valid_parent_indexes],
                           full_data$mass)
  full_data$color <- append(rep(0, length(scores)), full_data$color)
  list(full_data = full_data, parent_mass = parent_mass, index = index)
}




#' @export
#' @title Compute Molecular formula
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
compute_molecular_formulas4 <- function(mass_data, parent_ppm = 3,
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
  resultant_formulas <- DecomposeMasses(mass_data$ms2_matches$mz, mzs, parent_ppm, number_of_threads)
  mass_data$predicted_molecular_formulas <- resultant_formulas
  return(mass_data)

}





create_all_possible_formulas4 <- function(list_of_mz_int, parent_mass,
                                       parent_ppm, index) {
  parent_decomp <- decomposeMass(parent_mass, ppm = parent_ppm)
  valid_parent_indexes <- (which(parent_decomp$valid == "Valid"))
  invalid_indexes <- (which(parent_decomp$valid == "Invalid"))
  if (length(parent_decomp$formula) <= 0) {
    return(list(full_data = NA_character_, index = index))
  }
  if (length(valid_parent_indexes) == 1) {
    return(list(full_data = parent_decomp$formula[valid_parent_indexes],
                index = index))
  }
  if (length(parent_decomp$formula) == 1) {
    return(list(full_data = parent_decomp$formula[1], index = index))
  }
  # worse case scenario, use the invalid indexes to generate a value
  if (length(valid_parent_indexes) <= 0 && length(invalid_indexes) > 0) {
    valid_parent_indexes <- invalid_indexes
  }
  decomp_list <- vector("list", length(list_of_mz_int$mz))
  for (i in seq_along(list_of_mz_int$mz)) {
    if (parent_mass < list_of_mz_int$mz[[i]]) {
      next
    }
    decomp_list[[i]] <- decomposeIsotopes(list_of_mz_int$mz[[i]],
                                          list_of_mz_int$intensity[[i]])
  }
  full_data <- c()
  color_count <- 1
  for (i in seq_along(decomp_list)){
    if (is.null(decomp_list[[i]])) {
      next
    }
    valid_indexes <- (which(decomp_list[[i]]$valid == "Valid"))
    scores <- decomp_list[[i]]$score[valid_indexes]
    full_data$score <- append(full_data$score, scores)
    full_data$formula <- append(full_data$formula,
                                decomp_list[[i]]$formula[valid_indexes])
    full_data$mass <- append(full_data$mass,
                             decomp_list[[i]]$exactmass[valid_indexes])
    full_data$color <- c(full_data$color, rep(color_count, length(scores)))
    color_count <- color_count + 1
  }
  if (length(full_data$score) <= 0) {
    return(list(full_data = parent_decomp$formula[[1]],
                index = index))
  }
  scores <- parent_decomp$score[valid_parent_indexes]
  full_data$score <- append(scores, full_data$score)
  full_data$formula <- append(parent_decomp$formula[valid_parent_indexes],
                              full_data$formula)
  full_data$mass <- append(parent_decomp$exactmass[valid_parent_indexes],
                           full_data$mass)
  full_data$color <- append(rep(0, length(scores)), full_data$color)
  list(full_data = full_data, parent_mass = parent_mass, index = index)
}

#' @export
#' @title decomp
decomp_data <- function(mass, ppm = 3) {
  DecomposeMasses1(mass, ppm)
}




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
compute_molecular_formulas_other <- function(mass_data, parent_ppm = 3,
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
  resultant_formulas <- DecomposeMassesOther(mass_data$ms2_matches$mz, mzs, parent_ppm, number_of_threads)
  mass_data$predicted_molecular_formulas <- resultant_formulas
  return(mass_data)

}

