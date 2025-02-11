
#' @description
#' @export
diversity <- function(community_matrix, diversity_index) {
  if(!("matrix" %in% class(community_matrix))) {
    stop("Please ensure the community_matrix is a matrix")
  }
  # Checks for diversity_index on the c++ side
  CalculateDiversity(community_matrix, diversity_index)
}

#' @export
#' @description avgdist
averaged_subsampled_dissimilarity <- function(community_matrix, size, threshold, diversity_index = "bray", iterations = 1000) {
  return(FasterAvgDist(community_matrix, diversity_index, size, threshold, iterations))
}