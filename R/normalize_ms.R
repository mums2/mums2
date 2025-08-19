
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
  data
}
