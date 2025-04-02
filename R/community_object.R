#' @export
#' @title Create Community Matrix
#' @description
#' Takes the shared dataframe from clustur or a massdataset and converts it into a community matrix object
#' @param cluster_object the result of the `cluster_data()` function, or just a massdataset created from `convert_metaboscape2mass_dataset()` or `convert_mpactr_object_to_mass_data_set()`.
create_community_matrix_object <- function(cluster_object) {
  df <- get_abundance(cluster_object)
  samples <- unique(df$samples)
  combined_df <- data.frame(abund = df[which(df$samples == samples[[1]]), ]$abundance)

  for(i in 2:length(samples)) {
    combined_df <- cbind(combined_df, data.frame(abund = df[which(df$samples == samples[[i]]), ]$abundance))
  }

  combined_df <- t(as.matrix(combined_df))
  rownames(combined_df) <- samples
  obj <- CreateCommunityMatrix(combined_df)
  class(obj) <- c(class(obj), "community_object")
  return(obj)
}

create_community_matrix_object.mass_dataset <- function(data)
{
  data <- mass_data_set
  samples <- colnames(data@expression_data)
  matrix <- as.matrix(t(data@expression_data))
  rownames(matrix) <- samples
  mat <- CreateCommunityMatrix(matrix)
  class(mat) <- c(class(mat), "community_object")
  m <- get_community_matrix(mat)
}

create_community_matrix_object.clustur <- function(data)
{
  df <- get_abundance(data)
  samples <- unique(df$samples)
  combined_df <- data.frame(abund = df[which(df$samples == samples[[1]]), ]$abundance)

  for(i in 2:length(samples)) {
    combined_df <- cbind(combined_df, data.frame(abund = df[which(df$samples == samples[[i]]), ]$abundance))
  }

  combined_df <- t(as.matrix(combined_df))
  rownames(combined_df) <- samples
  obj <- CreateCommunityMatrix(combined_df)
  class(obj) <- c(class(obj), "community_object")
  return(obj)
}

#' @export
#' @title Get Community Matrix
#' @description
#' Returns the community `matrix` or the data that you used to create the object.
#' @param community_object the object created from the `create_community_object()` function.
get_community_matrix <- function(community_object) {
  if(!("community_object" %in% class(community_object))) {
    stop("Please ensure the community_object is created from the `create_community_object` function.")
  }
  return(GetCommunityMatrix(community_object))
}

#' @export
#' @title Print Community Object
#' @description
#' 3 function for print the community object
#' @param x the object created from the `create_community_object()` function.
#' @param ... other parameters that are included in the print function.
print.community_object <- function(x, ...) {
  print(get_community_matrix(x), ...)
}


