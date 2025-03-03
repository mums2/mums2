#' @title diversity_object
#' @export
#' @description diversity
diversity <- function(community_object, diversity_index){
  if(!("community_object" %in% class(community_object))) {
    stop("Please ensure the community_object is created from the `create_community_object` function.")
  }
  CalculateDiversityCommunityObject(community_object, diversity_index)
}

#' @title avgDist with Community object
#' @export
#' @description avgdist
averaged_subsampled_dissimilarity <- function(community_object, size, threshold, diversity_index = "bray", iterations = 1000) {
  return(FasterAvgDist(community_object, diversity_index, size, threshold, iterations))
}
