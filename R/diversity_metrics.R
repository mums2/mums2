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

#' @title Distance Shared
#' @export
#' @description dissimiliarty via beta diversity
#' @param community_object the object created from the `create_community_object()` function.
#' @param size the size you wish to rarefy your diversity matrix to.
#' @param threshold the threshold you want your species to reach before it is included
#' in the rarefaction sum.
#' @param diversity_index the diversity index you wish to calculate diversity, the two options are bray.
#' @param iterations the amount of times you wish to run your diversity metrics.
dist_shared <- function(community_object, size, threshold, diversity_index = "bray", iterations = 1000) {
  diversity_index_list <- c("bray", "jaccard", "soren", "hamming", "morisita", "thetayc")
  if(!("community_object" %in% class(community_object))) {
    stop("Please ensure the community_object is created from the `create_community_object` function.")
  }
  if(!(diversity_index %in% diversity_index_list)){
    stop(paste0("Please ensure your diversity index is one of the following values: ",
                paste(diversity_index_list, collapse = ', '))
    )
  }
  return(FasterAvgDist(community_object, diversity_index, size, threshold, iterations))
}


#' @title Alpha Diversity Summary
#' @export
#' @description alpha diversity
#' @param community_object the object created from the `create_community_object()` function.
#' @param size the size you wish to rarefy your diversity matrix to.
#' @param threshold the threshold you want your species to reach before it is included
#' in the rarefaction sum.
#' @param diversity_index the diversity index you wish to calculate diversity, the two options are
#' shannon or simpson.
#' @param iterations the amount of times you wish to run your diversity metrics.
alpha_summary <- function(community_object, size, threshold, diversity_index = "shannon", iterations = 1000) {
  diversity_index_list <- c("shannon", "simpson")
  if(!("community_object" %in% class(community_object))) {
    stop("Please ensure the community_object is created from the `create_community_object` function.")
  }

  if(!(diversity_index %in% diversity_index_list)){
    stop(paste0("Please ensure your diversity index is one of the following values: ",
                paste(diversity_index_list, collapse = ', '))
    )
  }
  return(FasterAvgDist(community_object, diversity_index, size, threshold, iterations))
}