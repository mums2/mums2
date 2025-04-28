#' @export
#' @title Convert distance data frame to a distance object
#' @description Converts the object generated from `dist_shared()` to a 
#' `dist` object.
#' @param distance_data_frame the object generated from the `dist_shared()` function.
community_object_to_distance_object <- function(distance_data_frame) {
  dist_object <- matrix(distance_data_frame$diversity,
         sqrt(nrow(distance_data_frame)), sqrt(nrow(distance_data_frame)))
  samples <- unique(c(distance_data_frame$firstSample, distance_data_frame$otherSample))
  rownames(dist_object) <- samples
  colnames(dist_object) <- samples
  return(as.dist(dist_object))
}