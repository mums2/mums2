#' @export
#' @title Cluster
#' @description
#' Clusters the data together
#' @param distance_df a distance df that was generated from the `distMs2()` function
#' @param mass_data_set your mass data set object
#' @param cluster_method a cluster method, there are five methods to choose from:
#' furthest, nearest, weighted, average, and opticlust. opticlust is the default
#' 
cluster_data <- function(distance_df, mass_data_set, cluster_method = "opticlust") {

  sparse_matrix <- create_sparse_matrix(distance_df$i, distance_df$j, distance_df$dist)

  # Create Count Table 
  count_table <- create_count_table(mass_data_set)

  # Create Distance Object
  dist <- read_dist(sparse_matrix, count_table, 0.3, F)

  # Cluster Data
  # results <- cluster(dist, 0.2, cluster_method)
  return(cluster(dist, 0.2, cluster_method))

}