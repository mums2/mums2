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
#' @param ms2_match_data your ms2-ms1 match data created from `ms2_ms1_compare()`
create_count_table <- function(ms2_match_data) {
  ms2_matches_compounds <- ms2_match_data$ms2_matches$ms1_compound_id
  peak_table <- ms2_match_data$ms1_data[ ,-c(2, 3, 4)]
  peak_table$cor <- NULL
  samples <- peak_table[which(peak_table$Compound %in% ms2_matches_compounds), ]
  return(data.frame(Representative_Sequence = samples$Compound, 
                    total = rowSums(samples[,-1]), samples[,-1], check.names = FALSE))
}