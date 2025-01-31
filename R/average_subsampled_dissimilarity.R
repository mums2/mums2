#' @export
#' @description avgdist
averaged_subsampled_dissimilarity <- function(community_matrix, size, threshold, diversity_index = "bray", iterations = 1000) {
  return(FasterAvgDist(community_matrix, diversity_index, size, threshold, iterations))
}