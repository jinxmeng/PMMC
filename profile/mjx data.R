# info Jinxin Meng, 20240301, 20240301 ----

setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/profile/")
pacman::p_load(dplyr, tidyr, tibble, purrr, ggplot2, ggpubr)

# data preprocess ----
source("F:/Code/R_func/transform_rc.R")
cvg <- read.delim("genomospecies.cvg", row.names = 1)
rc <- read.delim("genomospecies.rc", row.names = 1)
len <- read.delim("genomospecies.len", header = F, col.names = c("name", "len"))

## profile 按照宽度为50%过滤 ----
rc2 <- rc * map_df(cvg, \(x) ifelse(x < 0.5, 0, 1))
tpm <- rc2tpm(rc2[colSums(rc2)!=0], len)
saveRDS(rc, "genomospecies.rc.b50.rds")
saveRDS(tpm, "genomospecies.tpm.b50.rds")

tpm <- rc2tpm(rc, len)
saveRDS(tpm, "genomospecies.tpm.rds")

x = readRDS("genomospecies.tpm.b50.rds")

## 按照rc数过滤 ----
# rc2 <- rc * map_df(rc, \(x) ifelse(x < 100000, 0, 1))
# tpm <- rc2tpm(rc2[colSums(rc2)!=0], len)
# saveRDS(rc2, "genomospecies.rc.rds")
# saveRDS(tpm, "genomospecies.tpm.rds")


# variable ----
group <- read.delim("sample_group")

group_our <- group$group[group$group2 != ""] %>% unique
sample_our <- group %>% filter(group %in% group_our) %>% select(sample) %>% unlist(use.names = F)

group_region <- group %>% filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
  select(-sample) %>% unique.data.frame %>% 
  mutate(region = factor(region),
         region = forcats::fct_relevel(region, c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum"))) %>% 
  arrange(region, breed2) %>% select(group) %>% unlist(use.names = F)
sample_region <- group %>% filter(group %in% group_region) %>% select(sample) %>% unlist(use.names = F)

group_breed <- group$group[group$region == "Rectum"] %>% unique
sample_breed <- group %>% filter(group %in% group_breed) %>% select(sample) %>% unlist(use.names = F)

group_breed_pub <- group$group[group$group2 == ""] %>% unique
sample_breed_pub <- group %>% filter(group %in% group_breed_pub) %>% select(sample) %>% unlist(use.names = F)

group_breed_our <- group$group[group$group2 != "" & group$region == "Rectum"] %>% unique
sample_breed_our <- group %>% filter(group %in% group_breed_our) %>% select(sample) %>% unlist(use.names = F)

group_breed_bp <- group$group[group$class == "black" & group$region == "Rectum"] %>% unique
sample_breed_bp <- group %>% filter(group %in% group_breed_bp) %>% select(sample) %>% unlist(use.names = F)

region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")
save.image(file = "environ_info.RData")
