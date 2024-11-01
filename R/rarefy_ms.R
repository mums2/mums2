#' @importFrom tidyr uncount
#' @importFrom dplyr mutate
#' @importFrom dplyr "%>%"
#' @importFrom dplyr slice_head
#' @importFrom dplyr count
#' @importFrom dplyr filter
#' @export
rarefy_ms_tidy <- function(data, size, threshold) {
  random_pool <- data %>%
    tidyr::uncount(abund) %>%
    dplyr::mutate(mz = sample(mz))

  x <- size
  grand_total <- 0

  while (grand_total < size) {
    mz_counts <- random_pool %>%
      dplyr::slice_head(n = x) %>%
      dplyr::count(mz)
    grand_total <- sum(mz_counts$n[mz_counts$n > threshold])

    x <- x + 1
    # x <- x + (size - grand_total)
  }

  mz_counts %>%
    dplyr::filter(n > threshold)
}

#' @export
#' @importFrom tibble tibble
rarefy_ms_base <- function(data, size, threshold) {
  random_pool <- sample(rep(data$mz, data$abund))

  x <- size
  grand_total <- 0

  while(grand_total < size) {
    unfiltered <- table(random_pool[1:x])
    filtered <- unfiltered[unfiltered > threshold]
    grand_total <- sum(filtered)
    x <- x + (size - grand_total)
  }

  filtered_df <- tibble::tibble(mz = as.numeric(names(filtered)),
                                n = filtered)

  return(filtered_df)
}

#' @export
rarefy_ms_cpp <- function(data, size, threshold) {
  rarefyMs(data$mz, data$abund, size, threshold)
}

#' @export
rarefy_ms_cpp_2 <- function(data, size, threshold) {
  rarefyMs_2(data$mz, data$abund, size, threshold)
}
