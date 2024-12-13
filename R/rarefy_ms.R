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
#' @param data A `tibble` or `data.frame` where each row represets an MS1
#'  feature in a sample. The columns `mz` and `abund` are expected.
#' @param size The desired total sample intensity to subsample to.
#' @param threshold The individual feature threshold. Each subsampled feature
#'  must be >= this value to be retained.
#' 
#' @return A `tibble` or `data.frame` of rarefied feature intensities.
#' @export
#' 
#' @examples
#' set.seed(71)
#'
#' sample1 <- tibble::tibble(
#'    mz = seq(100, 1000, by = 100),
#'    abund = round(runif(10, 1000, 5e5)))
#' 
#' rarefy_ms(sample1, 25000, 1000)
#' 
rarefy_ms <- function(data, size, threshold) {
  rarefyMs(data$mz, data$abund, size, threshold)
}

#' Rarefy 
#'
#' @examples
rarefy_ms_generic <- function(data, size, threshold, feature_name = "mz") {
  rarefyMs_3(data[[feature_name]], data$abund, size, threshold)
}

#' Rarefy 
#'
#' @export
rarefy_four <- function(data, size, threshold, feature_name = "mz") {
  rarefyMs_4(data$mz, data$abund, size, threshold)
}

#' @export
benchmark <- function(){
 
  # microbenchmark::microbenchmark(CalculateAlphaDiversitySimpson(concentrated$mz, concentrated$abund, dilute_total, thresh, iterations = 1))
  # CalculateAlphaDiverstiyInt(concentrated$mz, concentrated$abund, dilute_total, thresh, iterations = 10),
  # times = 4)
}  

# microbenchmark::microbenchmark(CalculateAlphaDiversityShannon()(conc_rarefy$mz, conc_rarefy$abund, dilute_total, thresh, iterations = 1))
# CalculateBrayCurtisDissimilarity(list(dilute$mz, dilute$mz), list(dilute$abund, dilute$abund), dilute_total, thresh, iterations = 10)
