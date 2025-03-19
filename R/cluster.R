#' @export
#' @title Cluster
#' @description
#' Clusters the data together
#' 
cluster_data <- function(distance_df, peak_table, cluster_method = "opticlust") {

  sparse_matrix <- create_sparse_matrix(distance_df$i, distance_df$j, distance_df$dist)

  # Create Count Table 
  count_table <- create_count_table(peak_table)

  # Create Distance Object
  dist <- read_dist(sparse_matrix, count_table, 0.3, F)

  # Cluster Data
  return(cluster(dist, 0.2, cluster_method))

}