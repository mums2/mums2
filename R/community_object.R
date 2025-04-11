#' @export
#' @title Create Community Matrix
#' @description
#' Takes the shared dataframe from clustur or a massdataset and converts it into a community matrix object
#' @param data the result of the `cluster_data()` function, or just a mass_data object created from `ms2_ms1_compare()`.
create_community_matrix_object <- function(data) {
  return(UseMethod("create_community_matrix_object", data))
}

#' @export
#' @rdname create_community_matrix_object
create_community_matrix_object.mass_data <- function(data)
{
  data <- matched_data
  samples <- data$samples
  ms2_matches <-  data$ms2_matches$ms1_compound_id
  filtered_data <- data$ms1_data[which(data$ms1_data$Compound %in% ms2_matches), ][, samples, with = FALSE]
  matrix <- as.matrix(t(filtered_data))
  rownames(matrix) <- samples
  colnames(matrix) <- ms2_matches
  community_matrix <- CreateCommunityMatrix(matrix)
  class(community_matrix) <- c(class(community_matrix), "community_object")
  return(community_matrix)
}



#' @export
#' @rdname create_community_matrix_object
create_community_matrix_object.list <- function(data)
{
  data <- cluster_results
  df <- data$abundance
  samples <- unique(df$samples)
  combined_df <- data.frame(abund = df[which(df$samples == samples[[1]]), ]$abundance)

  for(i in 2:length(samples)) {
    combined_df <- cbind(combined_df, data.frame(abund = df[which(df$samples == samples[[i]]), ]$abundance))
  }

  combined_df <- t(as.matrix(combined_df))
  rownames(combined_df) <- samples
  colnames(combined_df) <- data$cluster$omu
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


