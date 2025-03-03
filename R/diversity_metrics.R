
#' @description
#' @export
diversity <- function(community_matrix, diversity_index) {
  if(!("matrix" %in% class(community_matrix))) {
    stop("Please ensure the community_matrix is a matrix")
  }
  # Checks for diversity_index on the c++ side
  CalculateDiversity(community_matrix, diversity_index)
}

#' @title diversity_object
#' @export
#' @description diversity
diversity_object <- function(community_object, diversity_index){
  CalculateDiversityCommunityObject(community_object, diversity_index)
}

#' @title avgDist
#' @export
#' @description avgdist
averaged_subsampled_dissimilarity <- function(community_matrix, size, threshold, diversity_index = "bray", iterations = 1000) {
  return(FasterAvgDist(community_matrix, diversity_index, size, threshold, iterations))
}

#' @title avgDist with Community object
#' @export
#' @description avgdist
averaged_subsampled_dissimilarity_object <- function(community_object, size, threshold, diversity_index = "bray", iterations = 1000) {
  return(FasterAvgDist2(community_object, diversity_index, size, threshold, iterations))
}

#' @title test1
#' @export
#' @description
#' A short description...
#' 
test_shuffle_No <- function(n)
{
  ShuffleVectorNoConversion(1:n)
}

#' @title test2
#' @export
#' @description
#' A short description...
#' 
test_shuffle <- function(n)
{
  ShuffleVectorConversion(1:n)
}

#' @title test2
#' @export
#' @description
#' A short description...
#' 
test_shuffle <- function(n)
{
 
}
