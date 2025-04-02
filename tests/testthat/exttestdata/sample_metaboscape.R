# This script generates a test feature table in Metaboscape format

# note: CSS is included, but as NA here to ensure this column is retained if included

library(tidyverse)

adducts <- c("ION=[M+H]+", "ION=[M+Na]+", "ION=[M+H+H2]3+", "ION=[M+K]+")

sample_ft <- data.frame("FEATURE_ID" = as.integer(1:10),
                        "RT" = sample(36:50, 10, replace = T),
                        "PEPMASS" = sample(400:1500, 10, replace = F),
                        "CSS" = rep(NA_real_, 10),
                        "ADDUCT" = sample(adducts, 10, replace = T),
                        "Blank_1_1" = c(rep(0, 4), round(runif(1, 100, 1500)), rep(0, 5)),
                        "Blank_1_2" = c(rep(0, 4), round(runif(1, 100, 1500)), rep(0, 5)),
                        "Control_1_1" = c(rep(0, 2), round(runif(6, 100, 1500)), rep(0, 2)),
                        "Control_1_2" = c(rep(0, 2), round(runif(6, 100, 1500)), rep(0, 2)),
                        "Exp_1_1" = c(rep(0, 2), round(runif(8, 100, 1500))),
                        "Exp_1_2" = c(rep(0, 2), round(runif(8, 100, 1500))))

write_csv(sample_ft, here::here("tests/testthat/exttestdata/sample_metaboscape.csv"))

sample_id <- sample_ft %>% 
  select(starts_with("Blank"), starts_with("Control"), starts_with("Exp")) %>% 
  colnames()

sample_info <- data.frame("sample_id" = sample_id,
                          "class" = c(rep("Blank", 2), rep("Sample", 4)),
                          "group" = c(rep("Solvent_blank", 2), rep("Control", 2), rep("Treatment1", 2)),
                          "injection.order" = sample(c(1:6), 6, replace = F))

write_csv(sample_info, here::here("tests/testthat/exttestdata/sample_metaboscape_metadata.csv"))
