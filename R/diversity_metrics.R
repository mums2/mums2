#' @title diversity_object
#' @export
#' @description diversity
#' @param community_object the object created from the `create_community_object()` function.
#' @param diversity_index the diversity index you wish to calculate diversity, the three options are
#' shannon, simpson, and bray.
diversity <- function(community_object, diversity_index){
  if(!("community_object" %in% class(community_object))) {
    stop("Please ensure the community_object is created from the `create_community_object` function.")
  }
  CalculateDiversityCommunityObject(community_object, diversity_index)
}

#' @title avgDist with Community object
#' @export
#' @description avgdist
#' @param community_object the object created from the `create_community_object()` function.
#' @param size the size you wish to rarefy your diversity matrix to.
#' @param threshold the threshold you want your species to reach before it is included
#' in the rarefaction sum.
#' @param diversity_index the diversity index you wish to calculate diversity, the three options are
#' shannon, simpson, and bray.
#' @param iterations the amount of times you wish to run your alpha or beta diversity metrics.
averaged_subsampled_dissimilarity <- function(community_object, size, threshold, diversity_index = "bray", iterations = 1000) {
  if(!("community_object" %in% class(community_object))) {
    stop("Please ensure the community_object is created from the `create_community_object` function.")
  }
  return(FasterAvgDist(community_object, diversity_index, size, threshold, iterations))
}

