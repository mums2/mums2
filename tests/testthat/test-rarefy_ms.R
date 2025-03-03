test_that("rarefy_ms returns the correct total", {
  # thresh <- 1000
  # set.seed(10)
  # concentrated <- tibble::tibble(
  #   mz = seq(100, 1000, by = 100),
  #   abund = round(runif(10, 1000, 5e5)))
  
  # dilute <- tibble::tibble(
  #   mz = concentrated$mz,
  #   abund = round(concentrated$abund / 100))
  
  # dilute_filter <- dplyr::filter(dilute, abund > thresh)

  # dilute_total <- sum(dilute_filter$abund)
  # # test <- as.matrix(t(concentrated))

  # # test <- test[2, ]
  # # test <- t(as.matrix(test))
  # # colnames(test) <- as.character(concentrated$mz)
  # #   rownames(test) <- "blank"
  # # div(test, "shannon")
  # # RarefactionCalculation(test, dilute_total, thresh)
  # conc_rarefy <- rarefy_ms(concentrated, dilute_total, thresh)
  # rare <- rarefy_four(concentrated, dilute_total, thresh)
  # # microbenchmark::microbenchmark(rarefy_ms(concentrated, dilute_total, thresh))
  # compare <- dplyr::full_join(concentrated, dilute_filter,
  #                      by = "mz", suffix = c(".conc", ".dil")) %>%
  #            dplyr::full_join(., conc_rarefy, by = "mz")
  
  # expect_equal(sum(compare$abund, na.rm = T), sum(compare$abund.dil, na.rm = T))
})
test_that("rarefy_ms returns the correct rowSum totals", {
  count_table <- test_path("exttestdata", "final.count_table")
  distances <- test_path("exttestdata", "final.dist")
  final_count <- read_count(count_table)
  final_dist <- read_dist(distances, final_count, 0.03, F)
  final_cluster <- cluster(final_dist, 0.03, "opticlust")
  size <- 400
  m <- create_community_matrix_object(final_cluster)
  resultant_matrix <- rarefy_ms(m, size, 10)
  expect_true(all(rowSums(resultant_matrix) >= size))

})
