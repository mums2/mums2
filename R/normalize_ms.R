
#' Normalize ms1 feature table
#'
#' @param .data A `mass_dataset` object with the slot `expression_data`
#' activated (see [massdataset::activate_mass_dataset()]). `expression_data`
#' is a `data.frame` with feature_id as rows and sample_id and columns.
#' Values are expected to represent absolute intensity for each feature in
#' the sample.
#' @param method A `character` to define which normalization method to apply.
#' Can be one of `c("rel_abund")`.
#'
#' @return A `data.frame` with feature_id as rows and sample_id and columns.
#' Values represent relative absolute intensity for each feature in the sample.
#' @export
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' dat <- tidyMassDemo
#'
#' # replace NA with 0 (or apply missing value imputation)
#' dat@expression_data <- dat@expression_data %>% replace(is.na(.), 0)
#'
#' dat_ms1_norm <- dat %>%
#'    massdataset::activate_mass_dataset("expression_data") %>%
#'    normalize_ms("rel_abund")
#'
#' # without normalization
#' dat@expression_data[1:5, 1:5]
#'
#' # with normalization
#' dat_ms1_norm@expression_data[1:5, 1:5]}
normalize_ms <- function(.data, method = NULL) {
  UseMethod("normalize_ms", .data)
}

#' @method normalize_ms mass_dataset
#' @importFrom cli cli_abort
#' @export
normalize_ms.mass_dataset <- function(.data, method = NULL) {

  if (length(.data@activated) == 0) {
    stop("activate your object using activate_mass_dataset first.\n")
  }

  # leave option to activate any slot for now, however this method
  # is specific to the expression_data slot of `mass_dataset` object
  x <- slot(object = .data, name = .data@activated)

  # check for na
  if (any(is.na(x) == TRUE)) {
    cli::cli_abort(c("Expression data contains NA values.
                      Clean data prior to between sample normalization."))
  }

  # Define renaming vector:
  methods <- c("rel_abund")

  # If method was not explicitly provided (is NULL), set to relative abundance
  if (is.null(method)){
    method <- methods[1]
  }

  if (method == "rel_abund") {
    expression_data_norm <- apply(x, 2, relative_abundance) %>%
      as.data.frame()
  }

  # expression_data_norm
  .data@expression_data <- expression_data_norm

  # update process info slot
  process_info <- .data@process_info

  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "mums2",
    function_name = "normalize_ms()",
    parameter = list(method = method),
    time = Sys.time()
  )

  if (all(names(process_info) != "normalize_ms")) {
    process_info$normalize_ms <- parameter
  } else {
    process_info$normalize_ms <- c(process_info$normalize_ms, parameter)
  }

  .data@process_info <- process_info

  return(.data)

}

################################
####        methods         ####
################################

#' Apply relative abundance normalization to ms1 feature table
#'
#' @details
#' `relative_abundance` normalizes feature intensities in a sample by
#' calculating the percent of each feature relative to the total
#' intensity of all features in the sample.
#'
#'
#' @param data A `data.frame` with feature_id as rows and sample_id and
#' columns. Values represent absolute intensity for each feature in the sample.
#'
#' @return A `data.frame` with feature_id as rows and sample_id and columns.
#' Values represent relative absolute intensity for each feature in the sample.
#'
#' @examples
#' sample_1 <- c(210.1333, 35.984, 21.264, 100.320, 3752.399)
#' sample_1_norm <- relative_abundance(sample_1)
#' @noRd 
relative_abundance <- function(data) {
  data <- data / sum(data)
  return(data)
}
