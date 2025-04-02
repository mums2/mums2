#' @export
#' @title create community matrix
#' @description
#' Takes the shared dataframe from clustur and converts it into a community matrix
#' @param cluster_object the result of the `cluster_data()` function. 
create_community_matrix <- function(cluster_object) {
  df <- get_abundance(cluster_object)
  samples <- unique(df$samples)
  combined_df <- data.frame(abund = df[which(df$samples == samples[[1]]), ]$abundance)

  for(i in 2:length(samples)) {
    combined_df <- cbind(combined_df, data.frame(abund = df[which(df$samples == samples[[i]]), ]$abundance))
  }

  combined_df <- t(as.matrix(combined_df))
  rownames(combined_df) <- samples
  return(combined_df)
}

#' @export
#' @title Create Count Table
#' @description
#' Creats a count table based with you mass data sets ms2 matches.
#' @param mass_data_set your mass data set object.
create_count_table <- function(mass_data_set) {
  ms2_matches <- mass_data_set@ms2_data[[1]]@variable_id
  samples <- mass_data_set@expression_data[which(rownames(mass_data_set@expression_data) %in% ms2_matches), ]
  return(data.frame(Representative_Sequence = rownames(samples), 
                    total = rowSums(samples), samples, check.names = FALSE))
  }

