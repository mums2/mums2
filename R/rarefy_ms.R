
#' @export
rarefy_ms <- function(data, size, threshold) {
  rarefyMs_2(data$mz, data$abund, size, threshold)
}
