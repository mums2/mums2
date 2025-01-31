#' Create `mass_dataset` from Metaboscape feature table
#'
#' @description
#' `convert_metaboscape2mass_dataset` takes a `data.frame` of sample
#' metadata (formatted with TidyMass requirements) and a `.csv` feature
#' table exported from Metaboscape to create the S4 `mass_dataset`
#'  object (Shen 2023).
#'
#' @details
#' This function allows users to import a Metaboscape feature table into
#'  TidyMass structure (`mass_dataset`).
#'
#' To create a `mass_dataset` object, we need at least three data frames:
#'  1. sample_info - a `data.frame` with, at minimum, the columns `sample_id`,
#'  `class`, `group`, and `injection.order`.
#'  1. variable_info - a `data.frame` with, at minimum, the columns
#'  `variable_id`, `mz`, and `rt`. The data frame holds our metadata for each
#'  MS1 feature in our feature table, where `variable_id` is the MS1 feature
#'  id, and `mz` and `rt` are the precursor mz and retention time, respectively,
#'  for each MS1 feature. Other feature information can be included, and is
#'  encouraged with the mums package function
#'  `convert_metaboscape2mass_dataset()` where any non-sample columns in the
#'  feature table .csv are retained in variable_info.
#'  1. expression_data - a `data.frame` where each row is an MS1 feature and
#'  each column is an injection sample.
#'
#' Information for `sample_info` of the `mass_dataset` object is provided in
#'  the `sample_info` argument. To format metadata for input as `sample_info`
#'  see https://massdataset.tidymass.org/articles/data_import_and_export.
#'
#' A Metaboscape feature table is expected for the argument `metaboscape_ft`,
#'  which  holds information for the `variable_info` and `expression_data`
#'  slots of the `mass_dataset` object. The function extracts all columns
#'  from the feature table that are not mean or absolute intensity values
#'  for groups/samples and stores them as `variable_info`. The
#'  `expression_data` is then extracted using sample_id from the `sample_info`
#'  `data.frame`. As a result, the sample_id in `sample_info` is expected
#'  to match sample column names in the Metaboscape feature table.
#'
#' @param metaboscape_ft a Metaboscape feature table as a `data.frame` or file
#' path.
#' @param sample_info a `data.frame` with, at minimum, the columns `sample_id`,
#' `class`, `group`, and `injection.order`.
#' @param sample_info_note optional, a `data.frame` with two columns `name` and
#' `meaning`. `name` should contain all column names in `sample_info` as
#' `character` `meaning` is a `character` vector describing the meaning for
#' the corresponding column.
#'
#' @return A `mass_dataset` object.
#'
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr select
#' @importFrom cli cli_abort
#' @importFrom massdataset create_mass_dataset
#'
#' @references
#' Shen X (2023). massdataset: massdataset. R package version 1.0.28,
#'  https://github.com/tidymass/massdataset.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # sample metadata
#' mjb_metadata
#' # Metaboscape feature table
#' mjb_metaboscape
#'
#' mjb_md <- convert_metaboscape2mass_dataset(mjb_metaboscape,
#'                                            mjb_metadata)}
#'
convert_metaboscape2mass_dataset <- function(metaboscape_ft,
                                             sample_info,
                                             sample_info_note = NULL) {
  cols <- c("sample_id", "class", "group", "injection.order")
  if (any(cols %in% colnames(sample_info) == FALSE)) {
    cli::cli_abort("{.cls {cols[which(!cols %in% colnames(sample_info))]}}
    are required columns that are not found in the provided
    sample_info. Please see function documentation for more
    details.")
  }

  if (!any(class(metaboscape_ft) == c("data.frame"))) {
    metaboscape_ft <- read.csv(metaboscape_ft)
  }

  if ("PEPMASS" %in% colnames(metaboscape_ft)) {
    if (!("ADDUCT" %in% colnames(metaboscape_ft))) {
      cli::cli_abort(c("Compund mass was found as PEPMASS, but column ADDUCT
      is not found in the data.frame. Mass cannot br converted
      to mz."))
    }
    cli::cli_alert(c("Compound mass was found as PEPMASS. This value will
    be converted to mz using ADDUCT"))
  }

  variable_info <- get_variable_info(metaboscape_ft, sample_info)

  expression_data <- metaboscape_ft %>%
    tibble::column_to_rownames(var = "FEATURE_ID") %>%
    dplyr::select(sample_info$sample_id)

  # check row names of expression_data match the variable_id column of
  #variable_info
  if (all(rownames(expression_data) != variable_info$variable_id)) {
    cli::cli_abort(c("feature ID (rownames) in expression data do not
    match variable_id in variable_info."))
  }

  # check the order of sample_id in sample_info and column names of expression
  #data are the same
  if (all(colnames(expression_data) != sample_info$sample_id)) {
    expression_data <- expression_data[, sample_info$sample_id]
  }

  if (is.null(sample_info_note)) {
    mass_dataset <- massdataset::create_mass_dataset(
                      expression_data = expression_data,
                      sample_info = sample_info,
                      variable_info = variable_info)
  } else {
    mass_dataset <- massdataset::create_mass_dataset(
                      expression_data = expression_data,
                      sample_info = sample_info,
                      variable_info = variable_info,
                      sample_info_note = sample_info_note)
  }

  return(mass_dataset)

}

#'
#' @importFrom dplyr select
#' @importFrom dplyr ends_with
#' @importFrom dplyr rename_all
#' @importFrom dplyr mutate
#' @importFrom dplyr case_when
#' @importFrom dplyr left_join
#'
#' @noRd
get_variable_info <- function(ft, sample_info) {
  file <- system.file("extdata/ion_masses", "DefinedIons.csv",
                      package = "mums2")
  ion_masses <- read.csv(file)

  ft %>%
    dplyr::select(-dplyr::ends_with("MeanIntensity"),
                  -sample_info$sample_id) %>%
    dplyr::rename_all(tolower) %>%
    dplyr::rename(variable_id = feature_id) %>%
    dplyr::mutate(variable_id = as.character(variable_id),
                  ion = gsub(".*\\[(.+)\\].*", "\\1", adduct),
                  charge_str = gsub(".*\\](.+)", "\\1", adduct),
                  charge = dplyr::case_when(charge_str == "+" ~ 1,
                                            charge_str == "2+" ~ 2,
                                            charge_str == "3+" ~ 3,
                                            TRUE ~ NA_real_)) %>%
    dplyr::left_join(ion_masses, dplyr::join_by("ion" == "IONS")) %>%
    dplyr::mutate(mz = (pepmass / charge) + MASS)
}


#' @export
convert_mpactr_object_to_mass_data_set <- function(mpactr_object) {
  dt <- as.data.table(get_peak_table(mpactr_object))
  meta <- get_meta_data(mpactr_object)
  # Get variable info
  variable_info <- dt[ , -c(unique(meta$Injection),"kmd"), with = FALSE]
  names(variable_info) <- c("variable_id", "mz", "rt")
   #sample_info
  sample_info <- data.frame(sample_id = meta$Injection, injection.order = 1:nrow(meta), class = meta$Sample_Code, 
                            group = meta$Biological_Group)
  
  # Get Expression data
  expression_data <- dt[ , unique(meta$Injection), with = FALSE]
  rownames(expression_data) <- variable_info$variable_id
  
  return(create_mass_dataset(expression_data = expression_data,
                      sample_info = sample_info,
                      variable_info = variable_info))
}

