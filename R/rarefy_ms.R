#' Rarefy MS1 Feature Table
#'
#' @description
#' `rarefy_ms()` performs a single subsampling of MS1 features in sample.
#'  Feature intensities are subsampled to the supplied `size` and accounts
#'  for intensity thresholds due to machine limits and background noise.
#'  Specifically, features whose abundance falls below the `threshold`
#'  after rarefying are removed. This allows for accurate representation
#'  of samples at different dilutions regardless of the desired
#'  submsampling `size`.
#'
#' @param community_object A `community_object`
#' @param size The desired total sample intensity to subsample to.
#' @param threshold The individual feature threshold. Each subsampled feature
#'  must be >= this value to be retained.
#' @param number_of_threads the number of threads you wish to use for this
#' calculation. Defaults to the number of threads on your computer.
#' @param seed the RNG (random number generator) seed you would like to use.
#' @return A `external_pointer` that references a community
#' matrix of rarefied feature intensities.
#' @export
#'
#' @examples
#' data <-
#'    import_all_data(peak_table =
#'                    mums2::mums2_example("botryllus_pt_small.csv"),
#'                    metadata =
#'                    mums2::mums2_example("boryillus_metadata.csv"),
#'                    format = "None")
#'
#'
#' matched_data <- ms2_ms1_compare(mums2_example("botryllus_v2.gnps.mgf"),
#'  data, 1, 6)
#'
#' community_object <- create_community_matrix_object(matched_data)
#' rarefy_ms(community_object, 400, 10)
#'
#' @return returns a `matrix` object that contains your rarefied data.
rarefy_ms <- function(community_object, size, threshold,
                      number_of_threads = detectCores(), seed = 123) {
  if (!inherits(community_object, "community_object")) {
    stop("Please ensure the community_object is created from the 
         `create_community_object` function.")
  }

  if (!is.numeric(size)) {
    stop("size must be numeric")
  }

  if (!is.numeric(threshold)) {
    stop("threshold must be numeric")
  }

  if (!is.numeric(number_of_threads)) {
    stop("number_of_threads must be numeric")
  }

  if (!is.numeric(seed)) {
    stop("seed must be numeric")
  }

  return(RarefactionCalculation(community_object, size,
                                threshold, number_of_threads, seed))
}
