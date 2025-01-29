# # dir <- "exttestdata"
# # file <- "demo_massdataset"
# # dat <- readRDS(test_path(dir, file))
# # data <- dat@ms2_data[[1]]@variable_id
# # dat_sub <- dat %>%
# #   massdataset::activate_mass_dataset("variable_info") %>%
# #   massdataset::filter(variable_id %in% data)
# # length(data)
# # dist_2 <- dist_ms2(dat_sub, 0.3, 2, gnps_params(0.5))
# # count_table <- data.frame(Representative_Sequences = data, total = rep(1, time = length(data)))
# # count_table_validated <- validate_count_table(count_table)
# # sparse_matrix <- clustur::create_sparse_matrix(dist_2$i, dist_2$j, dist_2$dist)
# # dist <- read_dist(sparse_matrix, count_table_validated, 0.3, F)
# # df <- get_distance_df(dist)
# # clust <- cluster(dist, 0.3, "opticlust", feature_column_name_to = "feature", bin_column_name_to = "omu")
# # r_file <- "database_data/PSU-MSMLS.msp"
# # psu_msmls <- massdatabase::read_msp_data(test_path(dir, r_file), source = "gnps")
# # cluster_list <- split_clusters_to_list(clust)
# # annotations <- annotate_ms2(dat_sub, psu_msmls,  gnps_params(0.5), 2, .7)
# # rarefy_data <- data.frame(mz = dat_sub@variable_info$mz)
# # rarefy_data_rounded <- data.frame(sapply(rarefy_data$mz, round, 2))
# # rarefy_data_corrected <- as.data.frame(table(rarefy_data_rounded))
# # names(rarefy_data_corrected) <- c("mz", "abund")
# # dat <- rarefy_ms(rarefy_data_corrected, sum(rarefy_data_corrected$abund), 5)
# data(BCI)
# # mean.avg.dist <- avgdist(dat, sample = 50, iterations = 10)
# # class(BCI)
# # which(cluster_list %in% annotations$query_ms1_id)
# # which(annotations$query_ms1_id %in% cluster_list[["omu438"]])
# # # Create Mass DataSet
# # sample_info_pos <- readr::read_csv("/Users/grejoh/Downloads/mpactR_test/dataset1_crc/metadata240418_allsamples_CompoundMeasurements.csv")
# # names(sample_info_pos)[1] <- "sample_id"
# # names(sample_info_pos)[3] <- "class"
# # variable_info_pos <- readr::read_csv("/Users/grejoh/Downloads/mpactR_test/dataset1_crc/peaktable240418_allsamples_CompoundMeasurements.csv", skip = 2)
# # expression_data_pos <- variable_info_pos[,4:ncol(variable_info_pos)]
# # variable_info_pos <- variable_info_pos[,1:3]
# # names(variable_info_pos) <- c("variable_id", "mz", "rt")
# # create_mass_dataset(
# #   expression_data = expression_data_pos,
# #   sample_info = sample_info_pos,
# #   variable_info = variable_info_pos
# # )

library(clustur)
# library(vegan)
library(labdsv)

#' @export
ta <- function()
{
  Test(conc_rarefy$abund)
}
#'
#' @export
rare_i <- function(data, size, threshold, feature_name = "mz") {
  rarefyMs_4(data$mz, data$abund, size, threshold)
}

# rrarefy(BCI, min(rowSums(BCI)))
# data(BCI)
# # Test the base functionality
# mean.avg.dist <- avgdist(BCI, sample = 50, iterations = 10)
# rare <- rrarefy(df_amazon, min(rowSums(df_amazon)))
v <- function(){
    data <- import_data(mpactr::example_path("cultures_metaboscape_peaktable.csv"), mpactr::example_path("cultures_metaboscape_metadata.csv"), "Metaboscape")
    data_filtered <- data |>
      filter_mispicked_ions(merge_peaks = TRUE, merge_method = "sum") |>
      filter_group(group_to_remove = "ANGDT") |>
      filter_cv(cv_threshold = 0.2, cv_param = "median")





  amazon_count <- read_count(example_path("amazon.sparse.count_table"))
  amazon_dist <- read_dist(example_path("amazon_column.dist"), amazon_count, 0.03, F)
  amazon_cluster <- cluster(amazon_dist, 0.03, "opticlust")
  samples <- unique(amazon_cluster$abundance$samples)
  amazon_sample_one <- amazon_cluster$abundance[which(amazon_cluster$abundance$samples == samples[[1]]), ]
  amazon_sample_two <- amazon_cluster$abundance[which(amazon_cluster$abundance$samples == samples[[2]]), ]
  combined_sample <- cbind(amazon_sample_one$abundance, amazon_sample_two$abundance)
  colnames(combined_sample) <- samples
  # Separate Samples and rename rows
  amazon_forest <- amazon_cluster$abundance[which(amazon_cluster$abundance$samples == "forest"), ]
  amazon_forest$bin <- 1:nrow(amazon_forest)
  amazon_pasture <- amazon_cluster$abundance[-which(amazon_cluster$abundance$samples == "forest"), ]
  amazon_pasture$bin <- 1:nrow(amazon_pasture)
  community_pasture <- matrify(amazon_pasture)
  community_forest <- matrify(amazon_forest)
  sum(amazon_pasture$abundance)
  community_mat <- matrify(amazon_cluster$abundance)

  amazon_shannon <- test_alpha_all(amazon_cluster)

  microbenchmark::microbenchmark(CalculateBrayCurtisDissimilarity(list(amazon_forest$bin, amazon_pasture$bin), 
    list(amazon_forest$abundance, amazon_pasture$abundance),30, 1, iterations=1000))
  microbenchmark::microbenchmark(CalculateAlphaDiversityShannon(amazon_forest$bin, amazon_forest$abundance, 30, 1, 1000))
  amazon_forest_generic <- amazon_forest
  colnames(amazon_forest_generic) <- c("sample" ,"mz", "abund")
  microbenchmark::microbenchmark(rarefy_ms(amazon_forest_generic, 30, 1))
  microbenchmark::microbenchmark(test_shannon())
  microbenchmark::microbenchmark(test_bray(1000))
}

get_sample <- function(sample, weight, sum){
  for(i in 1:sample){
    GetRandomNumberIndex(weight, length(weight), sum)
  }
}

final_dist_benchmark <- function(){
  final_count <- read_count("tests/testthat/exttestdata/final.count_table")
  final_dist <- read_dist("tests/testthat/exttestdata/final.dist", final_count, 0.03, F)
  final_cluster <- cluster(final_dist, 0.03, "opticlust")
  m <- create_community_matrix(final_cluster)
  rownames(m)
  diversity(m, "shannon")
  diversity(m, "simpson")
  diversity(m, "bray")
  r_m <- rarefaction(m, 400, 10)
  avg <- averaged_subsampled_dissimilarity(m, 400, 10)



  m <- prepare_for_rarefaction(final_cluster$abundance)
  microbenchmark::microbenchmark(rarefaction(m, 400, 10), rrarefy(m, 400))
  weight <- m[1, ]
  sum <- sum(weight)
  r_m <- rarefaction(m, 400, 5)
  GetRandomNumberIndex(weight, length(weight), sum)
  SomePaper(length(weight), 1, weight)
  avg_dist <- vegan::avgdist(m, 400)

  microbenchmark::microbenchmark(get_sample(400, weight, sum), SomePaper(length(weight), 400, weight))
  start_profiler("fast_avg_dist.out")
  my_avg_dist <- faster_avg_dist(m, "bray", 400, 10, 100)
  stop_profiler()
  avg_distance_function <- as.matrix(avgdist(m, 400, iterations = 100))
  comp_matrix <- abs(my_avg_dist - avg_distance_function)
  new_rarefaction(m, 400, 5)
  avg <- avgdist(m, 400, iterations = 1000)
  faster <- faster_avg_dist(m, "bray", 400, 10, 1000)
  microbenchmark::microbenchmark(averaged_subsampled_dissimilarity(m, 400, 10, "bray", 100), avgdist(m, 400, iterations = 100))
  sample_f3d2 <- final_cluster$abundance[which(final_cluster$abundance$samples == "F3D2"), ]
  sample_f3d2$bin <- 1:nrow(sample_f3d2)
  colnames(sample_f3d2)[2] <- "mz"
  colnames(sample_f3d2)[3] <- "abund"
  rarefy_ms(sample_f3d2, 10000, 100)
  community_mat <- matrify(sample_f3d2)
  sum(sample_f3d2$abund)
  diversity(final_cluster$abundance)
  microbenchmark::microbenchmark(
  CalculateAlphaDiversityShannon(sample_f3d2$mz, sample_f3d2$abund, 10000, 100), times = 10)
  set.seed(2)
  sample_f3d2_mat <- matrify(sample_f3d2)
  rarefy_ms(sample_f3d2, 10000, 100)
  microbenchmark::microbenchmark(test_shannon(sample_f3d2, 10000, 100), times = 5)
  microbenchmark::microbenchmark(diversity(rrarefy(sample_f3d2_mat, sample=10000)))
  microbenchmark::microbenchmark(
  test_alpha(sample_f3d2_mat), times = 10)
  diversity(rrarefy(community_pasture, sample=10000))
  d <- test_alpha_all(final_cluster)
  microbenchmark::microbenchmark(test_alpha_all(final_cluster), times = 5)
  rare <- test_new_rarefaction(final_cluster$abundance)
  CalculateDiversity(rare, "shannon")
  microbenchmark::microbenchmark(test_new_rarefaction(final_cluster$abundance), rrarefy(combined_df, 400))
  microbenchmark::microbenchmark(diversity(m, "shannon"), test_diversity_calc(m, "shannon"))
  start_profiler("rarefaction_test.out")
  microbenchmark::microbenchmark(new_rarefaction(m, 400, 5), rrarefy(m, 400))
  stop_profiler()
  m <- prepare_for_rarefaction(final_cluster$abundance)
  microbenchmark::microbenchmark(prepare_for_rarefaction(final_cluster$abundance))
  microbenchmark::microbenchmark(CalculateDiversity(m, "shannon"))
  microbenchmark::microbenchmark(diversity())
  microbenchmark::microbenchmark(RarefactionCalculation(m, 400, 10))
  p <- c(0.1,0.1,0.05,0.05, 0.2, 0.3, 0.05, 0.05, 0.05, 0.05)
  sum(p)
  length(p)
  unique(sample(1:10, size = 7, replace = T, prob = p))


  # Modifed Method
  faster <- faster_avg_dist(m, "bray", 400, 10, 1000)
  sample_names <- rownames(m)
  faster <- as.data.frame(faster)
  faster <- cbind(sample_names, faster)
  colnames(faster) <- c("names", sample_names)
  rownames(faster) <- sample_names
  fast <- pivot_longer(faster, -names)
  colnames(fast) <- c("rowname", "name", "value")

  # Vegan Method
  avg <- avgdist(m, 400, iterations = 1000)
  veg_df <- as.matrix(avg)
  samples <- colnames(veg_df)
  rownames(veg_df) <- samples
  veg_df <- cbind(samples, veg_df)
  colnames(veg_df) <- c("names", samples)
  veg <-  pivot_longer(as.data.frame(veg_df), -names)
  colnames(veg) <- c("rowname", "name", "value")
  veg$value <- as.double(veg$value)

  # Get mothur data
  mothur_matrix <- read.delim("tests/testthat/exttestdata/final.opti_mcc.braycurtis.0.03.square.ave.dist", 
  skip = 1,  col.names = letters[1:20], header = F)
  sample_names <- mothur_matrix[, 1]
  colnames(mothur_matrix) <- c("names" ,sample_names)
  mothur_1000 <- mothur_matrix %>% pivot_longer(-names)
  colnames(mothur_1000) <- c("rowname", "name", "value")


  inner_join(fast, mothur_1000, by = c("rowname", "name")) %>% ggplot(aes(x = value.x, y = value.y)) + 
    geom_point() + 
    geom_abline(intercept = 0, slope = 1) +
    ggtitle("Modified vs Mothur per 1000 iterations")

  inner_join(veg, mothur_1000, by = c("rowname", "name")) %>% ggplot(aes(x = value.x, y = value.y)) + 
    geom_point() + 
    geom_abline(intercept = 0, slope = 1) +
    ggtitle("Vegan vs Mothur per 1000 iterations")

  inner_join(fast, veg, by = c("rowname", "name")) %>% ggplot(aes(x = value.x, y = value.y)) + 
    geom_point() + 
    geom_abline(intercept = 0, slope = 1) +
    ggtitle("Modified vs Vegan per 1000 iterations")


  inner_join(veg, mothur_1000, by = c("rowname", "name")) %>% ggplot(aes(x = value.x, y = value.y)) + geom_point() + geom_abline(intercept = 0, slope = 1)


  mothur_matrix_longer <- mothur_matrix %>% pivot_longer(-names)
  inner_join(veg, mothur_matrix_longer, by = c("rowname", "name")) %>% ggplot(aes(x = value.x, y = value.y)) + geom_point() + geom_abline(intercept = 0, slope = 1)
  inner_join(veg, fast, by = c("rowname", "name"))

  inner_join(fast, mothur_matrix_longer, by = c("rowname", "name")) %>% ggplot(aes(x = value.x, y = value.y)) + geom_point() + geom_abline(intercept = 0, slope = 1)

  rownames_to_column(as.data.frame(faster)) %>% pivot_longer(-rowname) %>% filter(rowname < name) -> fast
  joined <- inner_join(veg, fast, by = c("rowname", "name"))
  joined <- inner_join(joined, mothur_1000, by = c("rowname", "name"))

}

rejection_sample <- function(n, s, p){
  k_i <- unique(sample(n, size = expected_items(n , s), T, p))
  l <- length(k_i)
  if(l >= s){
    return(k_i[[1:s]])               
  }

}


expected_items <- function(n,s)
{
  summation_to <- n
  summation_from <- n - s + 1
  sum <- 0
  for(i in summation_from:summation_to){
    sum <- sum + (1/i)
  }
  return(n * sum)
}

one_pass_sampling <- function(n, s, pi){
  # r  <- c()
  # for()
}
# microbenchmark::microbenchmark(wrswoR::sample_int_crank(1000, 200, 1:1000), new_rarefaction(m[1,] , 1000, 10))
# new_rarefaction(m[1,] , 10000, 100)

#' @export
new_rarefaction <- function(community_matrix, size, threshold){
  return(RarefactionCalculation(m, size, threshold))
}

#' @export
test_diversity_calc <- function(abundances, diversity_index) {
CalculateDiversity(abundances, diversity_index)
}

prepare_for_rarefaction <- function(df){
  samples <- unique(df$samples)
  combined_df <- data.frame(abund = df[which(df$samples == samples[[1]]), ]$abundance)

  for(i in 2:length(samples)) {
    combined_df <- cbind(combined_df, data.frame(abund = df[which(df$samples == samples[[i]]), ]$abundance))
  }

  combined_df <- t(as.matrix(combined_df))
  rownames(combined_df) <- samples
  return(combined_df)
}

test_new_rarefaction <- function(x)
{
  x <- final_cluster$abundance
  m <- prepare_for_rarefaction(x)
  samples <- unique(x$samples)
  rareified <- as.matrix(RarefactionCalculation(m, 400, 5))
  rareified <- t(rareified)
  row.names(rareified) <- samples
  return(as.matrix(rareified))
  
}

#'@export
faster_avg_dist <- function(community_matrix, diversity_index, sample, threshold, iterations)
{
  FasterAvgDist(community_matrix, diversity_index, sample, threshold, iterations)
}

test_alpha <- function(community_matrix) {
  sum <- 0
  sum <- sum + diversity(rrarefy(community_matrix, sample=10000))
  for(i in 1:1000) {
  }
  sum/1000
}
# matrix(1:3, 3, 3)

#' @export
benchmark_sample_without_replacement <- function()
{
  microbenchmark::microbenchmark(SomePaper(10, 1, p), GetRandomNumberIndex(p, 10, 1))
}

generate_sabund <- function(shared_df) {
  shared_df <- clustur::split_clusters_to_list(final_cluster)
  sabund <- lapply(shared_df, length)

}

#' @export
test_alpha_all <- function(shared_df) {
  # shared_df <- final_cluster
  samples <- unique(shared_df$abundance$samples)
  resultant_data_table <- data.frame(samples = samples)
  resultant_data_table$diversity <- 0
  for(i in 1:length(samples)){
    temp_df <- shared_df$abundance[which(shared_df$abundance$samples == samples[i]),]
    names(temp_df) <- c("samples", "mz", "abund")
    temp_df[[2]] <- 1:nrow(temp_df)
    resultant_data_table$diversity[[i]] <- CalculateAlphaDiversityShannon(temp_df$mz, temp_df$abund, sum(temp_df$abund)/2, 50)
  }
  return(resultant_data_table)
}
#' @export
test_shannon <- function(mat, size, threshold) {
  CalculateAlphaDiversityShannon(mat$mz, mat$abund, size, threshold)
}

#' @export
test_bray <- function(iters = 10)
{
  CalculateBrayCurtisDissimilarity(list(amazon_forest$bin, amazon_pasture$bin), 
  list(amazon_forest$abundance, amazon_pasture$abundance),30, 1, iters)
}

test_diversity <- function(x){
  browser()
  diversity(x)
}
test_rrarefy <- function(x) {
  browser()
  rrarefy(x, 100)
}
# microbenchmark::microbenchmark(rrarefy(community_mat, 5))
# microbenchmark::microbenchmark(diversity(community_mat, index = "shannon"))
# microbenchmark::microbenchmark(avgdist(community_mat, 30, iterations = 1000))
# test_bray(iters = 1000)
# microbenchmark::microbenchmark(CalculateBrayCurtisDissimilarity(list(amazon_forest$bin, amazon_pasture$bin), 
#   list(amazon_forest$abundance, amazon_pasture$abundance),30, 1, iterations=1000))
# # diversity(df_amazon) 
# df_amazon <- reshape2::dcast(amazon_cluster$abundance, samples ~ bin)
# df_amazon <- matrify(amazon_cluster$abundance)
# avg_dist <- avgdist(df_amazon, sample = 50, iterations = 1000)
# diversity(df_amazon, index = "simpson")
# df <- reshape2::dcast(amazon_cluster$abundance, samples ~ bin)
# modified_abundance_df <- amazon_cluster$abundance
# modified_abundance_df$omu <- 1:length(modified_abundance_df$bin)
# mod_df_rarefy <-  modified_abundance_df[]
# colnames(modified_abundance_df)[2] <- "mz"
# colnames(modified_abundance_df)[3] <- "abund"
# sum_abund <- sum(modified_abundance_df$abund)
# dat_2 <- rarefy_ms(modified_abundance_df[,2:3], round(sum_abund/2), 4)

# df["samples"] <- 1
# avgdist(df_amazon, sample = 50, iterations = 10, binary = T)
# dat <- rarefy_ms(rarefy_data_corrected, sum(rarefy_data_corrected$abund), 5)
# microbenchmark::microbenchmark(GetRandomNumberIndex(numbers_sample, 100), RcppSample(numbers_sample, number_weights), times = 1000)
# set.seed(40)


# func <- function(){
#   data(BCI)

#   browser()
#   bray <- vegdist(BCI, "bray")
#   rrarefy(BCI, 50)
#   mean.avg.dist <- avgdist(BCI, sample = 50, iterations = 100)
# }
# vegdist(sub_amazon, "bray")
# bray <- vegdist(BCI, "bray", binary = TRUE)
# func()

# #' @export 
# func <- function() {
#   microbenchmark::microbenchmark(CalculateAlphaDiversityInt(concentrated$mz, concentrated$abund, dilute_total, thresh, iterations = 1000), times = 5)
# }

#' @export
cpp <- function(iter = 1) {
  CalculateAlphaDiversityShannon(concentrated$mz, concentrated$abund, dilute_total, thresh, iterations = iter)
}


# cpp()

# cpp()

# CalculateAlphaDiversityInt(concentrated$mz, concentrated$abund, dilute_total, thresh, iterations = 1000)
# microbenchmark::microbenchmark(cpp)
# # microbenchmark::microbenchmark(CalculateAlphaDiversityShannon(concentrated$mz, concentrated$abund, dilute_total, thresh, iterations = 1), times = 10)
# microbenchmark::microbenchmark(diversity(rrarefy(m, 100), "shannon"))

#' @export
cpp <- function(iter = 1) {
  CalculateAlphaDiversityShannon(concentrated$mz, concentrated$abund, dilute_total, thresh, iterations = iter)
}

#' @export
bray <- function(iter = 1) {
  CalculateBrayCurtisDissimilarity(list(concentrated$mz), list(concentrated$abund), dilute_total, thresh, iterations = iter)
}



# # cpp()
# microbenchmark::microbenchmark(cpp(1000), times = 10)
fun <- function() {
  concentrated$samples <- rep("no_group", times = 10)
  test <- data.frame(sample = concentrated$samples, mz = concentrated$mz, abund = concentrated$abund)  
  m <- matrify(test)
  # microbenchmark::microbenchmark(rrarefy(m, sample = 25011))
  microbenchmark::microbenchmark(rrarefy(m, sample = dilute_total),  rarefy_four(concentrated, dilute_total, thresh))
}

f <- function()
{
  browser()
  diversity(m, "simpson")
}
# sub.sample
# subsample.h/subsample.cpp
# summary.shared command
# 

# diversity(m, "shannon")
# f()

# benchmark()
# conc_two <- concentrated
# conc_two$mz <- as.character(conc_two$mz)
# conc_rarefy <- rarefy_ms_generic(conc_two, dilute_total, thresh)
# name <- "mz"
# conc_two[[name]]

# microbenchmark::microbenchmark(CalculateAlphaDiverstiy(conc_two$mz, conc_two$abund, dilute_total, thresh, iterations = 10),
# CalculateAlphaDiverstiyInt(concentrated$mz, concentrated$abund, dilute_total, thresh, iterations = 1000),
# times = 5)

# func <- benchmark2() {
#   microbenchmark::microbenchmark(CalculateAlphaDiversityInt(concentrated$mz, concentrated$abund, dilute_total, thresh, iterations = 1000), times = 5)
# }

# microbenchmark::microbenchmark(CalculateAlphaDiversityInt(concentrated$mz, concentrated$abund, dilute_total, thresh, iterations = 1000), times = 5)
# # microbenchmark::microbenchmark(vegan::rarefy())

# data <- import_data(example("cultures_peak_table.csv"),
#   example("cultures_metadata.csv"),
#   format = "Progenesis"
# )


# data_mpactr <- filter_mispicked_ions(data,
#   ringwin = 0.5,
#   isowin = 0.01,
#   trwin = 0.005,
#   max_iso_shift = 3,
#   merge_peaks = TRUE,
#   merge_method = "sum"
# )
# data_mpactr <- filter_group(data_mpactr, 0.01, "Solvent_Blank", TRUE)

# data_mpactr_copy <- filter_insource_ions(data_mpactr,
#                                          cluster_threshold = 0.95,
#                                          copy_object = TRUE)

# data_mpactr <- filter_insource_ions(data_mpactr, cluster_threshold = 0.95)

# print(data_mpactr)
# filtered <- get_peak_table(data_mpactr)[1:5, 1:5]
# pt <- get_peak_table(data_mpactr) %>%
#   select(Compound, starts_with("102423"), starts_with("102623")) %>%
#   pivot_longer(-Compound, names_to = "sample") %>% #10152
#   mutate(value = as.integer(value)) %>%
#   mutate(total = sum(value), .by = c("sample")) %>%
#   filter(total > 0) %>% #8,883
#   mutate(total = sum(value), .by = c("Compound")) %>% #8,883
#   filter(total != 0) %>% #8,883
#   select(-total) %>%
#   filter(!str_detect(sample, "Media")) %>%
#   mutate(Compound = str_c("ion", Compound, sep = "_"))
# head(pt)
# pt_df_temp <- as.data.frame(pt)
# mat_pt <- matrify(pt_df_temp)


# pt_df <- pt %>%
#   pivot_wider(names_from = "Compound", values_from = "value", values_fill = 0)
#   column_to_rownames(var = "sample")

#   dim(pt_df)

# rare <- rarefy(mat_pt, 284)
# microbenchmark::microbenchmark(rarefy(mat_pt, 284), times = 5)
# res <- avgdist(mat_pt, 10, iterations = 1000)
# pt_df_2 <- as.data.frame(pt)
# mat_pt <- labdsv::matrify(pt_df_2)
# print(mat_pt)

# We cant really avg dist with mz unless we add the column, do we need the mz? What if we cluster the data, what then?
# TODO: Create a bray curtis calculator
# TODO: Keep researching about bias in alpha diversity
# TODO: Find ways to replicate the data in vegan (alpha/beta diversity. We have the formula)


# test_log <- function(){
#   product <- 1
#   for(i in 1:100){
#     product <- product * i
#   }
#   return(log(product))
# }

# test_log_2 <- function(){
#   sum <- 0
#   for(i in 1:100){
#     sum <- sum + log(i)
#   }
#   return(sum)
# }
# test_log()
# test_log_2()
# microbenchmark::microbenchmark(test_log(), test_log_2(), times = 10)

# sum(rarefy_four(concentrated, dilute_total, thresh)$abund)
# ls <- vector("numeric", 1000)
# ls_data <- list()

# for(i in 1:1000){
#   rarefy_four(concentrated, dilute_total, thresh)
#   # summation <- sum(dat$abund) 
#   # # if(summation < dilute_total) {
#   # #   ls_data <- c(ls_data, dat)
#   # #   next
#   # # }
#   # ls[i] <- summation
# }
# # sum(ls_data$abund)
# # max(ls)
# # min(ls)

# 1 - dilute_total/max(ls)
# 1 - (min(ls))/dilute_total




# TODO Make rarefaction and diversity functions faster than mothur and vegan
# TODO Finish the pipeline for mums2
# TODO Benchmark package


# SomePaper(100, 10, p)
