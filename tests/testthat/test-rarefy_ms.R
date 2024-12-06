test_that("rarefy_ms returns the correct total", {
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
  microbenchmark::microbenchmark(rarefy_ms(concentrated, dilute_total, thresh))
  compare <- dplyr::full_join(concentrated, dilute_filter,
                       by = "mz", suffix = c(".conc", ".dil")) %>%
             dplyr::full_join(., conc_rarefy, by = "mz")
  
  expect_equal(sum(compare$abund, na.rm = T), sum(compare$abund.dil, na.rm = T))
})


# benchmark()
# conc_two <- concentrated
# conc_two$mz <- as.character(conc_two$mz)
# conc_rarefy <- rarefy_ms_generic(conc_two, dilute_total, thresh)
# name <- "mz"
# conc_two[[name]]

# microbenchmark::microbenchmark(CalculateAlphaDiverstiy(conc_two$mz, conc_two$abund, dilute_total, thresh, iterations = 10),
# CalculateAlphaDiverstiyInt(concentrated$mz, concentrated$abund, dilute_total, thresh, iterations = 10),
# times = 10)

# # microbenchmark::microbenchmark(vegan::rarefy())

data <- import_data(example("cultures_peak_table.csv"),
  example("cultures_metadata.csv"),
  format = "Progenesis"
)


data_mpactr <- filter_mispicked_ions(data,
  ringwin = 0.5,
  isowin = 0.01,
  trwin = 0.005,
  max_iso_shift = 3,
  merge_peaks = TRUE,
  merge_method = "sum"
)
data_mpactr <- filter_group(data_mpactr, 0.01, "Solvent_Blank", TRUE)

data_mpactr_copy <- filter_insource_ions(data_mpactr,
                                         cluster_threshold = 0.95,
                                         copy_object = TRUE)

data_mpactr <- filter_insource_ions(data_mpactr, cluster_threshold = 0.95)

print(data_mpactr)
filtered <- get_peak_table(data_mpactr)[1:5, 1:5]
pt <- get_peak_table(data_mpactr) %>%
  select(Compound, starts_with("102423"), starts_with("102623")) %>%
  pivot_longer(-Compound, names_to = "sample") %>% #10152
  mutate(value = as.integer(value)) %>%
  mutate(total = sum(value), .by = c("sample")) %>%
  filter(total > 0) %>% #8,883
  mutate(total = sum(value), .by = c("Compound")) %>% #8,883
  filter(total != 0) %>% #8,883
  select(-total) %>%
  filter(!str_detect(sample, "Media")) %>%
  mutate(Compound = str_c("ion", Compound, sep = "_"))
head(pt)
pt_df_temp <- as.data.frame(pt)
mat_pt <- matrify(pt_df_temp)


pt_df <- pt %>%
  pivot_wider(names_from = "Compound", values_from = "value", values_fill = 0)
  column_to_rownames(var = "sample")

  dim(pt_df)

rare <- rarefy(mat_pt, 284)
microbenchmark::microbenchmark(rarefy(mat_pt, 284), times = 5)
res <- avgdist(mat_pt, 10, iterations = 1000)
pt_df_2 <- as.data.frame(pt)
mat_pt <- labdsv::matrify(pt_df_2)
print(mat_pt)

# We cant really avg dist with mz unless we add the column, do we need the mz? What if we cluster the data, what then?
# TODO: Create a bray curtis calculator
# TODO: Keep researching about bias in alpha diversity
# TODO: Find ways to replicate the data in vegan (alpha/beta diversity. We have the formula)