#' Rarefy MS1 Feature Table
#' 
#' @description
#' `rarefy_ms()` performs a single subsampling of MS1 features in sample.
#'  Feature intensities are subsampled to the supplied `size` and accounts
#'  for intesnity thresholds due to machine limits and background noise.
#'  Specifically, features whose abundance falls below the `threshold`
#'  after rarefying are removed. This allows for accurate represetation
#'  of samples at diffrent dilutions regardless of the desired
#'  submsampling `size`.
#' 
#' @param community_object A `community_object`
#' @param size The desired total sample intensity to subsample to.
#' @param threshold The individual feature threshold. Each subsampled feature
#'  must be >= this value to be retained.
#' 
#' @return A `matrix` of rarefied feature intensities.
#' @export
#' 
#' @examples
#' set.seed(71)
#'
#' sample1 <- tibble::tibble(
#'    mz = seq(100, 1000, by = 100),
#'    abund = round(runif(10, 1000, 5e5)))
#' 
#'
#' 
rarefy_ms <- function(community_object, size, threshold) {
  if(!("community_object" %in% class(community_object))) {
    stop("Please ensure the community_object is created from the `create_community_object` function.")
  }
  return(RarefactionCalculation(community_object, size, threshold))
}
