test_that("rarefy_ms_tidy returns the correct total", {
  set.seed(3) # failed on seed 103
  thresh <- 1000

  concentrated <- tibble::tibble(
    mz = seq(100, 1000, by = 100),
    abund = round(runif(10, 1000, 5e5)))
  
  dilute <- tibble::tibble(
    mz = concentrated$mz,
    abund = round(concentrated$abund / 100))
  
  dilute_filter <- dplyr::filter(dilute, abund > thresh)

  dilute_total <- sum(dilute_filter$abund)

  conc_rarefy <- rarefy_ms_tidy(concentrated, dilute_total, thresh)

  compare <- dplyr::full_join(concentrated, dilute_filter,
                       by = "mz", suffix = c(".conc", ".dil")) %>%
             dplyr::full_join(., conc_rarefy, by = "mz")
  
  expect_equal(sum(compare$n, na.rm = T), sum(compare$abund.dil, na.rm = T))
})

test_that("rarefy_ms_base returns the correct total", {
  set.seed(6)
  thresh <- 1000

  concentrated <- tibble::tibble(
    mz = seq(100, 1000, by = 100),
    abund = round(runif(10, 1000, 5e5)))
  
  dilute <- tibble::tibble(
    mz = concentrated$mz,
    abund = round(concentrated$abund / 100))
  
  dilute_filter <- dplyr::filter(dilute, abund > thresh)

  dilute_total <- sum(dilute_filter$abund)

  conc_rarefy <- rarefy_ms_base(concentrated, dilute_total, thresh)

  compare <- dplyr::full_join(concentrated, dilute_filter,
                       by = "mz", suffix = c(".conc", ".dil")) %>%
             dplyr::full_join(., conc_rarefy, by = "mz")
  
  expect_equal(sum(compare$n, na.rm = T), sum(compare$abund.dil, na.rm = T))
})

test_that("rarefy_ms_cpp returns the correct total", {
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

  # conc_rarefy <- rarefyMs_2(concentrated$mz, concentrated$abund, dilute_total, thresh)
  conc_rarefy <- rarefy_ms_cpp(concentrated, dilute_total, thresh)

  compare <- dplyr::full_join(concentrated, dilute_filter,
                       by = "mz", suffix = c(".conc", ".dil")) %>%
             dplyr::full_join(., conc_rarefy, by = "mz")
  
  expect_equal(sum(compare$abund, na.rm = T), sum(compare$abund.dil, na.rm = T))
})

test_that("rarefy_ms_cpp_2 returns the correct total", {
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
  
  conc_rarefy <- rarefy_ms_cpp_2(concentrated, dilute_total, thresh)

  compare <- dplyr::full_join(concentrated, dilute_filter,
                       by = "mz", suffix = c(".conc", ".dil")) %>%
             dplyr::full_join(., conc_rarefy, by = "mz")
  
  expect_equal(sum(compare$abund, na.rm = T), sum(compare$abund.dil, na.rm = T))
})


# microbenchmark::microbenchmark(rarefy_ms_cpp(concentrated, dilute_total, thresh), 
# rarefy_ms_cpp_2(concentrated, dilute_total, thresh), times=10)
# ##### benchmarking
# load_all()
# thresh <- 1000
# concentrated <- tibble::tibble(
#   mz = seq(100, 1000, by = 100),
#   abund = round(runif(10, 1000, 5e5)))

# dilute <- tibble::tibble(
#   mz = concentrated$mz,
#   abund = round(concentrated$abund / 100))

# dilute_filter <- dplyr::filter(dilute, abund > thresh)

# dilute_total <- sum(dilute_filter$abund)

# library("mums2")
# microbenchmark::microbenchmark(rarefy_ms_tidy(concentrated, dilute_total, thresh),
# rarefy_ms_base(concentrated, dilute_total, thresh),
# rarefy_ms_cpp(concentrated, dilute_total, thresh),
# rarefy_ms_cpp_2(concentrated, dilute_total, thresh),
#   times = 100
# )
#
# Note: to accurately assess c++ functions, you need to install() and not load_all()
# as load_all() includes debugging information that signigicantly impacts speed.
#
# Following install(), the c++ functions are faster than the r functions (3 - 191 times faster).
# The base R function outperforms the tidy R function by 3917 milisections, on average.
# rarefy_ms_cpp_2 is faster than rarefy_ms_cpp by 11 miliseconds.
#
# Performace order from fastest to slowest is now: cpp2 < cpp < base < tidy.

