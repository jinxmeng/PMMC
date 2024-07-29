# Jinxin Meng, 20240305, 20240714 ------------------

setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/compos/")
pacman::p_load(dplyr, tidyr, tibble, purrr, ggplot2, ggpubr)
library(randomForest)

# rf -----------------------
source("/code/R_func/randomforest.R")
source("/code/R_func/taxa.R")
source("F:/Code/R_func/difference_analysis.R")

group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/genomospecies.tpm.b50.rds")
taxonomy <- read.delim("../profile/genome_taxonomy.txt") %>% 
  mutate(phylum = gsub("_\\w$", "", phylum),
         family = gsub("_\\w$", "", family),
         genus = gsub("_\\w$", "", genus))

region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")
breed_order <- c("JL-BSWP","JL-LXP","JL-SLBP","JL-DLLW","HB-MSP","HL-MP","HL2-MP","HN-YNBP",
                 "IM-IMBP","LN-BJBP","LN-JSBP","LN-LDBP","LN-THNP","LN2-LDBP","SD-LWP","SX-LLBP")

group_x <- group %>%
  filter(group2 != "" & region == "Rectum") %>% 
  select(sample, group = breed2)

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

# family ----------------

data <- taxa_trans(profile_x, taxonomy, to = "genus", out_all = T, transRA = T)
breeds <- setdiff(breed_order, "JL-DLLW")

out <- list()
for (i in breeds) {
  groups <- c(i, "JL-DLLW")
  group_i <- filter(group_x, group %in% groups)
  profile_i <- select(data, any_of(group_i$sample)) %>% 
    filter(rowSums(.) != 0)
  rf <- rf_var_rank(profile_i, group_i, seed = 2024, ntree = 1000)
  out[[i]] <- rf
}

map(out, \(x) x %>% filter(MeanDecreaseAccuracy > 0) %>% pull(1)) %>%
  unlist() %>% 
  table %>% sort()

group_i <- group_x %>% 
  mutate(group = ifelse(group == "JL-DLLW", "A", "B"))
profile_i <- select(data, any_of(group_i$sample)) %>% 
  filter(rowSums(.) != 0)
rf <- rf_var_rank(profile_i, group_i, seed = 2024, ntree = 1000)


