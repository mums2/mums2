#' @export
#' @description
#' Takes the shared dataframe from clustur and converts it into a community matrix
#' 
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
#' Creats a count table based on your peak table
#' 
create_count_table <- function(peak_table) {
  sample_cols <- colnames(peak_table)[5:ncol(peak_table)]
    
  count_table <- data.frame(Representative_Sequence = 1:nrow(peak_table))
  count_table$sum <- rowSums(peak_table[, .SD, .SDcols = sample_cols])
  count_table <- cbind(count_table, peak_table[, .SD, .SDcols = sample_cols])
  colnames(count_table)[2] <- "total"
  count_table$Representative_Sequence <- peak_table$Compound
  count_table$Representative_Sequence <- as.character(count_table$Representative_Sequence)
  return(count_table)
 }
