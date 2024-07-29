# Jinxin Meng, 20240305, 20240701 ------------------

setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/compos/")
pacman::p_load(dplyr, tidyr, tibble, purrr, ggplot2, ggpubr)

# prok composition phylum -----------------

source("/code/R_func/taxa.R")
source("/code/R_func/calcu_diff.R")

group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/genomospecies.tpm.b50.rds")
taxonomy <- read.delim("../profile/genome_taxonomy.txt") %>% 
  mutate(phylum = gsub("_\\w$", "", phylum),
         family = gsub("_\\w$", "", family))

region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")
breed_order <- c("JL-BSWP","JL-LXP","JL-SLBP","JL-DLLW","HB-MSP","HL-MP","HL2-MP","HN-YNBP",
                 "IM-IMBP","LN-BJBP","LN-JSBP","LN-LDBP","LN-THNP","LN2-LDBP","SD-LWP","SX-LLBP")

group_x <- group %>%
  filter(group2 != "") %>% 
  select(sample, group, breed = breed2, region) %>% 
  mutate(region = factor(region, region_order),
         breed = factor(breed, breed_order))

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

# sample level
data <- taxa_trans(profile_x, taxonomy, group_x, to = "phylum", top_n = 10,
                   transRA = T, other_name = "p__Other")

plot_data <- data %>% 
  rownames_to_column(var = "taxa") %>% 
  gather(key = "sample", value = "value", -taxa) %>% 
  left_join(group_x, by = "sample")

group_order <- plot_data %>%
  arrange(region, breed) %>% 
  pull(group) %>% 
  as.character() %>% 
  unique()

sample_order <- plot_data %>% 
  filter(taxa == "p__Bacteroidota") %>% 
  mutate(group = factor(group, group_order)) %>% 
  arrange(group, desc(value)) %>% 
  pull(sample)

plot_data <- mutate(plot_data, sample = factor(sample, sample_order))

colors <- c("#bebada","#ffffb3","#8dd3c7","#FF9D9A","#fee08b",
            "#a7c9e0","#de77ae","#fabfd2","#b8e186","#d9d9d9")

ggplot(plot_data, aes(x = sample, y = value, fill = taxa)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 2),
        legend.text = element_text(size = 8, color = "#000000", face = "italic"),
        panel.grid = element_blank())

ggsave("compos/prok.compos.phylum.by_smps.pdf", height = 5, width = 15)

# group level
data <- taxa_trans(profile_x, taxonomy, group_x, to = "phylum", out_all = T, 
                   transRA = T, smp2grp = T, method = "mean", other_name = "p__Other")

rownames_to_column(data, var = "name") %>% 
  write.table("compos/prok.compos.phylum.by_grps.tsv", sep = "\t", quote = F, row.names = F)

data <- taxa_trans(profile_x, taxonomy, group_x, to = "phylum", top_n = 10,
                   transRA = T, smp2grp = T, method = "mean",  other_name = "p__Other")

plot_data <- data %>% 
  rownames_to_column(var = "taxa") %>% 
  gather(key = "group", value = "value", -taxa) %>% 
  left_join(group_x %>% select(-sample) %>% unique, by = "group") %>% 
  mutate(group = factor(group, group_order))

ggplot(plot_data, aes(x = group, y = value, fill = taxa)) +
  geom_bar(stat = "identity", color = "black", width = .7, size = .3) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
        legend.text = element_text(size = 8, color = "#000000", face = "italic"),
        panel.grid = element_blank())

ggsave("compos/prok.compos.phylum.by_grps.pdf", height = 5, width = 10)

# prok composition family -----------------

# sample level
data <- taxa_trans(profile_x, taxonomy, group_x, to = "family", top_n = 10,
                   transRA = T, other_name = "f__Other")

plot_data <- data %>% 
  rownames_to_column(var = "taxa") %>% 
  gather(key = "sample", value = "value", -taxa) %>% 
  left_join(group_x, by = "sample")

group_order <- plot_data %>%
  arrange(region, breed) %>% 
  pull(group) %>% 
  as.character() %>% 
  unique()

sample_order <- plot_data %>% 
  filter(taxa == "f__Lachnospiraceae") %>% 
  mutate(group = factor(group, group_order)) %>% 
  arrange(group, desc(value)) %>% 
  pull(sample)

plot_data <- mutate(plot_data, sample = factor(sample, sample_order))

colors <- c("#bebada","#de77ae","#8dd3c7","#fee08b","#FF9D9A",
            "#a7c9e0","#ffffb3","#fabfd2","#b8e186","#d9d9d9")

ggplot(plot_data, aes(x = sample, y = value, fill = taxa)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 2),
        legend.text = element_text(size = 8, color = "#000000", face = "italic"),
        panel.grid = element_blank())

ggsave("compos/prok.compos.family.by_smps.pdf", height = 5, width = 15)

# group level
data <- taxa_trans(profile_x, taxonomy, group_x, to = "family", out_all = T, 
                   transRA = T, smp2grp = T, method = "mean", other_name = "f__Other")

rownames_to_column(data, var = "name") %>% 
  write.table("compos/prok.compos.family.by_grps.tsv", sep = "\t", quote = F, row.names = F)

data <- taxa_trans(profile_x, taxonomy, group_x, to = "family", top_n = 10,
                   transRA = T, smp2grp = T, method = "mean", other_name = "f__Other")

plot_data <- data %>% 
  rownames_to_column(var = "taxa") %>% 
  gather(key = "group", value = "value", -taxa) %>% 
  left_join(group_x %>% select(-sample) %>% unique, by = "group") %>% 
  mutate(group = factor(group, group_order))

ggplot(plot_data, aes(x = group, y = value, fill = taxa)) +
  geom_bar(stat = "identity", color = "black", width = .7, size = .3) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
        legend.text = element_text(size = 8, color = "#000000", face = "italic"),
        panel.grid = element_blank())

ggsave("compos/prok.compos.phylum.by_grps.pdf", height = 5, width = 10)

# phage composition family ------------------------------------------------

source("/code/R_func/taxa.R")
source("/code/R_func/calcu_diff.R")

group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/final_vOTUs.tpm.rds")
metadata <- read.delim("../phage_genome/data.phage_genome_metadata.txt")
taxonomy <- select(metadata, name, family, type)

region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")
breed_order <- c("JL-BSWP","JL-LXP","JL-SLBP","JL-DLLW","HB-MSP","HL-MP","HL2-MP","HN-YNBP",
                 "IM-IMBP","LN-BJBP","LN-JSBP","LN-LDBP","LN-THNP","LN2-LDBP","SD-LWP","SX-LLBP")

group_x <- group %>%
  filter(group2 != "") %>% 
  select(sample, group, breed = breed2, region) %>% 
  mutate(region = factor(region, region_order),
         breed = factor(breed, breed_order))

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

# sample level
data <- taxa_trans(profile_x, taxonomy, group_x, to = "family", top_n = 10,
                   transRA = T, other_name = "Other")

plot_data <- data %>% 
  rownames_to_column(var = "taxa") %>% 
  gather(key = "sample", value = "value", -taxa) %>% 
  left_join(group_x, by = "sample")

group_order <- plot_data %>%
  arrange(region, breed) %>% 
  pull(group) %>% 
  as.character() %>% 
  unique()

sample_order <- plot_data %>% 
  filter(taxa == "Siphoviridae") %>% 
  mutate(group = factor(group, group_order)) %>% 
  arrange(group, desc(value)) %>% 
  pull(sample)

plot_data <- mutate(plot_data, sample = factor(sample, sample_order))

colors <- c("#bebada","#de77ae","#8dd3c7","#fee08b","#FF9D9A",
            "#a7c9e0","#ffffb3","#fabfd2","#b8e186","#d9d9d9")
colors <- c("#bebada","#fee08b","#8dd3c7","#a7c9e0","#de77ae",
            "#fabfd2","#b8e186","#FF9D9A","#ffffb3","#d9d9d9")

ggplot(plot_data, aes(x = sample, y = value, fill = taxa)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 2),
        legend.text = element_text(size = 8, color = "#000000", face = "italic"),
        panel.grid = element_blank())

ggsave("compos/phage.compos.family.by_smps.pdf", height = 5, width = 15)

# group level
data <- taxa_trans(profile_x, taxonomy, group_x, to = "family", out_all = T, 
                   transRA = T, smp2grp = T, method = "mean", other_name = "Other")

rownames_to_column(data, var = "name") %>% 
  write.table("compos/phage.compos.family.by_grps.tsv", sep = "\t", quote = F, row.names = F)

data <- taxa_trans(profile_x, taxonomy, group_x, to = "family", top_n = 10,
                   transRA = T, smp2grp = T, method = "mean", other_name = "Other")

plot_data <- data %>% 
  rownames_to_column(var = "taxa") %>% 
  gather(key = "group", value = "value", -taxa) %>% 
  left_join(group_x %>% select(-sample) %>% unique, by = "group") %>% 
  mutate(group = factor(group, group_order))

ggplot(plot_data, aes(x = group, y = value, fill = taxa)) +
  geom_bar(stat = "identity", color = "black", width = .7, size = .3) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
        legend.text = element_text(size = 8, color = "#000000", face = "italic"),
        panel.grid = element_blank())

ggsave("compos/phage.compos.phylum.by_grps.pdf", height = 5, width = 10)
