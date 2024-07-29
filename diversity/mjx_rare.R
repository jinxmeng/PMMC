# Jinxin Meng, 20240301, 20240701 --------------------------

setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/diversity/")
pacman::p_load(dplyr, tidyr, tibble, purrr, ggplot2, ggpubr, ggpmisc)

# prok ------------------------------------------------------

source("/code/R_func/rare_curve.R")
group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/genomospecies.tpm.b50.rds")

# our data
group_x <- filter(group, group2 != "")

profile_x <- select(profile, any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

data <- calcu_specaccum(profile_x, permutations = 99)
write.table(data, "rare/prok.rare.specaccum.476_smps.tsv", sep = "\t", row.names = F, quote = F)

ggline(data, x = "name", y = "obs", point.size = NA, 
       xlab = "Number of samples", ylab = "Number of genomospecies")
ggsave("rare/prok.rare.specaccum.476_smps.pdf", width = 5, height = 4)

# data <- calcu_specaccum_by_group(profile_x, group_x)
# write.table(data, "prok.rare.specaccum.476_smps.by_group.tsv", sep = "\t", row.names = F, quote = F)

# ggline(data, x = "name", y = "obs", group = "group", point.size = NA, 
#        xlab = "Number of samples", ylab = "Number of genomospecies")
# ggsave("prok.rare.specaccum.476_smps.by_group.pdf", width = 5, height = 4)

# public data
group_x <- filter(group, group2 == "") 

profile_x <- select(profile, any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

data <- calcu_specaccum(profile_x, permutations = 99)
write.table(data, "rare/prok.rare.specaccum.pub_smps.tsv", sep = "\t", row.names = F, quote = F)

ggline(data, x = "name", y = "obs", point.size = NA, 
       xlab = "Number of samples", ylab = "Number of genomospecies")
ggsave("rare/prok.rare.specaccum.pub_smps.pdf", width = 5, height = 4)

# data <- calcu_specaccum_by_group(profile_x, group_x)
# write.table(data, "prok.rare.specaccum.pub_smps.by_group.tsv", sep = "\t", row.names = F, quote = F)

# ggline(data, x = "name", y = "obs", group = "group", point.size = NA, 
#        xlab = "Number of samples", ylab = "Number of genomospecies")
# ggsave("prok.rare.specaccum.pub_smps.by_group.pdf", width = 5, height = 4)

# phage ----------------------------------------------- 

profile <- readRDS("../profile/final_vOTUs.tpm.rds")

# our data
group_x <- filter(group, group2 != "")

profile_x <- select(profile, any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

data <- calcu_specaccum(profile_x, permutations = 99)
write.table(data, "rare/phage.rare.specaccum.476_smps.tsv", sep = "\t", row.names = F, quote = F)

ggline(data, x = "name", y = "obs", point.size = NA, 
       xlab = "Number of samples", ylab = "Number of vOTUs")
ggsave("rare/phage.rare.specaccum.476_smps.pdf", width = 5, height = 4)

# data <- calcu_specaccum_by_group(profile_x, group_x)
# write.table(data, "phage.rare.specaccum.476_smps.by_group.tsv", sep = "\t", row.names = F, quote = F)
# 
# ggline(data, x = "name", y = "obs", group = "group", point.size = NA, 
#        xlab = "Number of samples", ylab = "Number of vOTUs")
# ggsave("phage.rare.specaccum.476_smps.by_group.pdf", width = 5, height = 4)

# public data
group_x <- filter(group, group2 == "")

profile_x <- select(profile, any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

data <- calcu_specaccum(profile_x, permutations = 99)
write.table(data, "rare/phage.rare.specaccum.pub_smps.tsv", sep = "\t", row.names = F, quote = F)

ggline(data, x = "name", y = "obs", point.size = NA, 
       xlab = "Number of samples", ylab = "Number of vOTUs")
ggsave("rare/phage.rare.specaccum.pub_smps.pdf", width = 5, height = 4)

# data <- calcu_specaccum_by_group(profile_x, group_x)
# write.table(data, "phage.rare.specaccum.pub_smps.by_group.tsv", sep = "\t", row.names = F, quote = F)
# 
# ggline(data, x = "name", y = "obs", group = "group", point.size = NA, 
#        xlab = "Number of samples", ylab = "Number of vOTUs")
# ggsave("phage.rare.specaccum.pub_smps.by_group.pdf", width = 5, height = 4)
