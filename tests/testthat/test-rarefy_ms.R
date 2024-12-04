test_that("rarefy_ms returns the correct total", {
  set.seed(71) # failed on seed 92
  thresh <- 1000

  concentrated <- tibble::tibble(
    mz = seq(100, 1000, by = 100),
    abund = round(runif(10, 1000, 5e5)))
  
  dilute <- tibble::tibble(
    mz = concentrated$mz,
    abund = round(concentrated$abund / 100))
  
  dilute_filter <- dplyr::filter(dilute, abund > thresh)

  dilute_total <- sum(dilute_filter$abund)

  conc_rarefy <- rarefy_ms(concentrated, dilute_total, thresh)

  compare <- dplyr::full_join(concentrated, dilute_filter,
                       by = "mz", suffix = c(".conc", ".dil")) %>%
             dplyr::full_join(., conc_rarefy, by = "mz")
  
  expect_equal(sum(compare$abund, na.rm = T), sum(compare$abund.dil, na.rm = T))
})
