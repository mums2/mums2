dir <- "exttestdata"
file <- "demo_massdataset"
dat <- readRDS(test_path(dir, file))
data <- dat@ms2_data[[1]]@variable_id
dat_sub <- dat %>%
  massdataset::activate_mass_dataset("variable_info") %>%
  massdataset::filter(variable_id %in% data)
length(data)
dist_2 <- dist_ms2(dat_sub, 0.3, 2, gnps_params(0.5))
count_table <- data.frame(Representative_Sequences = data, total = rep(1, time = length(data)))
count_table_validated <- validate_count_table(count_table)
sparse_matrix <- clustur::create_sparse_matrix(dist_2$i, dist_2$j, dist_2$dist)
dist <- read_dist(sparse_matrix, count_table_validated, 0.3, F)
df <- get_distance_df(dist)
clust <- cluster(dist, 0.3, "opticlust", feature_column_name_to = "feature", bin_column_name_to = "omu")
r_file <- "database_data/PSU-MSMLS.msp"
psu_msmls <- massdatabase::read_msp_data(test_path(dir, r_file), source = "gnps")
cluster_list <- split_clusters_to_list(clust)
annotations <- annotate_ms2(dat_sub, psu_msmls,  gnps_params(0.5), 2, .7)
rarefy_data <- data.frame(mz = dat_sub@variable_info$mz)
rarefy_data_rounded <- data.frame(sapply(rarefy_data$mz, round, 2))
rarefy_data_corrected <- as.data.frame(table(rarefy_data_rounded))
names(rarefy_data_corrected) <- c("mz", "abund")
dat <- rarefy_ms(rarefy_data_corrected, sum(rarefy_data_corrected$abund), 5)
data(BCI)
mean.avg.dist <- avgdist(dat, sample = 50, iterations = 10)
class(BCI)
which(cluster_list %in% annotations$query_ms1_id)
which(annotations$query_ms1_id %in% cluster_list[["omu438"]])
# Create Mass DataSet
sample_info_pos <- readr::read_csv("/Users/grejoh/Downloads/mpactR_test/dataset1_crc/metadata240418_allsamples_CompoundMeasurements.csv")
names(sample_info_pos)[1] <- "sample_id"
names(sample_info_pos)[3] <- "class"
variable_info_pos <- readr::read_csv("/Users/grejoh/Downloads/mpactR_test/dataset1_crc/peaktable240418_allsamples_CompoundMeasurements.csv", skip = 2)
expression_data_pos <- variable_info_pos[,4:ncol(variable_info_pos)]
variable_info_pos <- variable_info_pos[,1:3]
names(variable_info_pos) <- c("variable_id", "mz", "rt")
create_mass_dataset(
  expression_data = expression_data_pos,
  sample_info = sample_info_pos,
  variable_info = variable_info_pos
)


amazon_count <- read_count(example_path("amazon.sparse.count_table"))
amazon_dist <- read_dist(example_path("amazon_column.dist"), amazon_count, 0.03, F)
amazon_cluster <- cluster(amazon_dist, 0.03, "opticlust")
diversity(df_amazon)

df_amazon <- reshape2::dcast(amazon_cluster$abundance, samples ~ bin)
diversity(df_amazon)
df <- reshape2::dcast(clust$abundance, samples ~ omu)
modified_abundance_df <- clust$abundance
modified_abundance_df$omu <- 1:length(modified_abundance_df$omu)
colnames(modified_abundance_df)[2] <- "mz"
colnames(modified_abundance_df)[3] <- "abund"
sum_abund <- sum(modified_abundance_df$abund)
dat_2 <- rarefy_ms(modified_abundance_df[,2:3], round(sum_abund/2), 4)

df["samples"] <- 1
avgdist(df, sample = 1, iterations = 10)
dat <- rarefy_ms(rarefy_data_corrected, sum(rarefy_data_corrected$abund), 5)
# microbenchmark::microbenchmark(GetRandomNumberIndex(numbers_sample, 100), RcppSample(numbers_sample, number_weights), times = 1000)
# RcppSample(numbers_sample, number_weights, 3)
# set.seed(40)
# GetRandomNumberIndex(c(10,20,30,40), 100)