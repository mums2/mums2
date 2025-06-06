#' @export
#' @title Cluster
#' @description
#' Clusters the data together
#' @param distance_df a distance df that was generated from the `distMs2()` function.
#' @param ms2_match_data your mass data set object generated from `ms2_ms1_compare()`.
#' @param cutoff the cutoff value you wish to cluster to.
#' @param cluster_method a cluster method, there are five methods to choose from:
#' furthest, nearest, weighted, average, and opticlust.
#' @return a shared `data.frame` displaying all the clustered and abundance data.
cluster_data <- function(distance_df, ms2_match_data, cutoff = 0.3, cluster_method = "opticlust") {

  sparse_matrix <- create_sparse_matrix(distance_df$i, distance_df$j, distance_df$dist)

  # Create Count Table 
  count_table <- create_count_table(ms2_match_data)

  # Create Distance Object
  dist <- read_dist(sparse_matrix, count_table, cutoff, FALSE)

  # Cluster Data
  return(cluster(dist, cutoff, cluster_method, bin_column_name_to = "omu"))
}