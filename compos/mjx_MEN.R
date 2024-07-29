# Jinxin Meng, 20240305, 20240712 ------------------

setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/compos/")
pacman::p_load(dplyr, tidyr, tibble, purrr, ggplot2, ggpubr)

# prok MEN for_breed family ------------------------

source("F:/Code/R_func/taxa.R")
source("F:/Code/R_func/calcu_MEN.R")

group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/genomospecies.tpm.b50.rds")
taxonomy <- read.delim("../profile/genome_taxonomy.txt") %>% 
  mutate(phylum = gsub("_\\w$", "", phylum),
         family = gsub("_\\w$", "", family))

group_x <- group %>%
  filter(group2 != "" & region == "Rectum") %>% 
  select(sample, group = breed2, class)

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0) %>% 
  taxa_trans(taxonomy, to = "family", out_all = T) %>% 
  filter(rownames(.) != "f__Unknown")

breeds <- unique(group_x$group)
plot_list <- list()
attr_list <- list()
for (i in breeds) {
  profile_i <- select(profile_x, any_of(group_x$sample[group_x$group == i]))
  profile_i <- profile_i[rowSums(profile_i) != 0, colSums(profile_i) != 0]
  metadata <- data.frame(family = rownames(profile_i)) %>% 
    left_join(select(taxonomy, phylum, family) %>% unique(), by = "family") %>% 
    column_to_rownames("family")
  
  data <- calcu_MEN(profile_i, metadata, title = i) %>% suppressWarnings()
  plot_list[[i]] <- data$gplot
  attr_list[[i]] <- data$network_attr
}

# 网络图
cowplot::plot_grid(plotlist = plot_list, nrow = 3)
ggsave("MEN/prok.network.for_breed.family.pdf", width = 30, height = 15)

# 网络属性
attr_data <- attr_list %>% 
  map2(., names(.), \(x, y) rename(x, !!y := value)) %>% 
  reduce(\(x, y) merge(x, y, by = "name"))
write.table(attr_data, "MEN/prok.network.for_breed.family.attr.tsv", sep = "\t", row.names = F, quote = F)

plot_data <- gather(attr_data, key = "group", value = "value", -name) %>% 
  left_join(select(group_x, group, class) %>% unique, by = "group")

ggbarplot(plot_data, x = "group", y = "value", fill = "group", legend = "none",
          position = position_dodge(), x.text.angle = 90, size = .4, ggtheme = theme_bw(),
          font.tickslab = c(8, "plain", "#000000"), xlab = "", ylab = "") +
  geom_hline(aes(yintercept = value),
             filter(plot_data, group == "JL-DLLW"),
             lty = "dashed", linewidth = .4) +
  facet_wrap(~factor(name), scales = "free_y") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", color = "#000000"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = .4, color = "#000000"),
        panel.background = element_rect(linewidth = .4, color = "#000000"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1/2)
ggsave("MEN/prok.network.for_breed.family.attr.barplot.pdf", width = 15, height = 9)

# phage MEN for_breed family ------------------------

source("F:/Code/R_func/taxa.R")
source("F:/Code/R_func/calcu_MEN.R")

group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/final_vOTUs.tpm.rds")
taxonomy <- read.delim("../phage_genome/data.phage_genome_metadata.txt") %>% 
  select(name, family)

group_x <- group %>%
  filter(group2 != "" & region == "Rectum") %>% 
  select(sample, group = breed2, class)

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0) %>% 
  taxa_trans(taxonomy, to = "family", out_all = T)

breeds <- unique(group_x$group)
plot_list <- list()
attr_list <- list()
for (i in breeds) {
  profile_i <- select(profile_x, any_of(group_x$sample[group_x$group == i]))
  profile_i <- profile_i[rowSums(profile_i) != 0, colSums(profile_i) != 0]
  metadata <- data.frame(family = rownames(profile_i)) %>% 
    mutate(phylum = paste0("phy_", 1:nrow(.))) %>% 
    column_to_rownames("family")
  
  data <- calcu_MEN(profile_i, metadata, title = i) %>% suppressWarnings()
  plot_list[[i]] <- data$gplot
  attr_list[[i]] <- data$network_attr
}

# 网络图
cowplot::plot_grid(plotlist = plot_list, nrow = 3)
ggsave("MEN/phage.network.for_breed.family.pdf", width = 30, height = 15)

# 网络属性
attr_data <- attr_list %>% 
  map2(., names(.), \(x, y) rename(x, !!y := value)) %>% 
  reduce(\(x, y) merge(x, y, by = "name"))
write.table(attr_data, "MEN/phage.network.for_breed.family.attr.tsv", sep = "\t", row.names = F, quote = F)

plot_data <- gather(attr_data, key = "group", value = "value", -name) %>% 
  left_join(select(group_x, group, class) %>% unique, by = "group") %>% 
  mutate(value = ifelse(grepl("[A-Za-z]", value), 0, value)) %>% 
  mutate(value = round(as.numeric(value), 2))

ggbarplot(plot_data, x = "group", y = "value", fill = "group", legend = "none",
          position = position_dodge(), x.text.angle = 90, size = .4, ggtheme = theme_bw(),
          font.tickslab = c(8, "plain", "#000000"), xlab = "", ylab = "") +
  geom_hline(aes(yintercept = value),
             filter(plot_data, group == "JL-DLLW"),
             lty = "dashed", linewidth = .4) +
  facet_wrap(~factor(name), scales = "free_y") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", color = "#000000"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = .4, color = "#000000"),
        panel.background = element_rect(linewidth = .4, color = "#000000"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1/2)
ggsave("MEN/phage.network.for_breed.family.attr.barplot.pdf", width = 15, height = 9)

# combine MEN for_breed family ------------------------

source("F:/Code/R_func/taxa.R")
source("F:/Code/R_func/calcu_MEN.R")

group <- read.delim("../profile/sample_group")
profile_x <- readRDS("../profile/genomospecies.tpm.b50.rds")
taxonomy_x <- read.delim("../profile/genome_taxonomy.txt") %>% 
  mutate(phylum = gsub("_\\w$", "", phylum),
         family = gsub("_\\w$", "", family))
profile_y <- readRDS("../profile/final_vOTUs.tpm.rds")
taxonomy_y <- read.delim("../phage_genome/data.phage_genome_metadata.txt") %>% 
  select(name, family)

group_x <- group %>%
  filter(group2 != "" & region == "Rectum") %>% 
  select(sample, group = breed2, class)

profile_x <- rbind(
  profile_x %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0) %>% 
  taxa_trans(taxonomy_x, to = "family", out_all = T),
  profile_y %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0) %>% 
  taxa_trans(taxonomy_y, to = "family", out_all = T)
  )

breeds <- unique(group_x$group)
plot_list <- list()
attr_list <- list()
for (i in breeds) {
  profile_i <- select(profile_x, any_of(group_x$sample[group_x$group == i]))
  profile_i <- profile_i[rowSums(profile_i) != 0, colSums(profile_i) != 0]
  metadata <- data.frame(family = rownames(profile_i)) %>% 
    mutate(phylum = paste0("phy_", 1:nrow(.))) %>% 
    column_to_rownames("family")
  
  data <- calcu_MEN(profile_i, metadata, title = i) %>% suppressWarnings()
  plot_list[[i]] <- data$gplot
  attr_list[[i]] <- data$network_attr
}

# 网络图
cowplot::plot_grid(plotlist = plot_list, nrow = 3)
ggsave("MEN/combine.network.for_breed.family.pdf", width = 30, height = 15)

# 网络属性
attr_data <- attr_list %>% 
  map2(., names(.), \(x, y) rename(x, !!y := value)) %>% 
  reduce(\(x, y) merge(x, y, by = "name"))
write.table(attr_data, "MEN/combine.network.for_breed.family.attr.tsv", sep = "\t", row.names = F, quote = F)

plot_data <- gather(attr_data, key = "group", value = "value", -name) %>% 
  left_join(select(group_x, group, class) %>% unique, by = "group") %>% 
  mutate(value = ifelse(grepl("[A-Za-z]", value), 0, value)) %>% 
  mutate(value = round(as.numeric(value), 2))

ggbarplot(plot_data, x = "group", y = "value", fill = "group", legend = "none",
          position = position_dodge(), x.text.angle = 90, size = .4, ggtheme = theme_bw(),
          font.tickslab = c(8, "plain", "#000000"), xlab = "", ylab = "") +
  geom_hline(aes(yintercept = value),
             filter(plot_data, group == "JL-DLLW"),
             lty = "dashed", linewidth = .4) +
  facet_wrap(~factor(name), scales = "free_y") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", color = "#000000"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = .4, color = "#000000"),
        panel.background = element_rect(linewidth = .4, color = "#000000"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1/2)
ggsave("MEN/combine.network.for_breed.family.attr.barplot.pdf", width = 15, height = 9)

# prok MEN for_breed genomospecies ------------------------

source("F:/Code/R_func/taxa.R")
source("F:/Code/R_func/calcu_MEN.R")

group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/genomospecies.tpm.b50.rds")

group_x <- group %>%
  filter(group2 != "" & region == "Rectum") %>% 
  select(sample, group = breed2, class)

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

# breeds <- unique(group_x$group)
# plot_list <- list()
# attr_list <- list()
# for (i in breeds) {
#   profile_i <- select(profile_x, any_of(group_x$sample[group_x$group == i]))
#   profile_i <- profile_i[rowSums(profile_i) != 0, colSums(profile_i) != 0]
#   metadata <- data.frame(taxa = rownames(profile_i)) %>% 
#     mutate(phylum = paste0("phy_", 1:nrow(.))) %>% 
#     column_to_rownames("taxa")
#   
#   data <- calcu_MEN(profile_i, metadata, title = i) %>% suppressWarnings()
#   plot_list[[i]] <- data$gplot
#   attr_list[[i]] <- data$network_attr
# }
# save(plot_list, attr_list, file = "prok.for_breed.genomospecies.RData")
load("MEN/prok.for_breed.genomospecies.RData")

# 网络图
cowplot::plot_grid(plotlist = plot_list, nrow = 3)
ggsave("MEN/prok.network.for_breed.genomospecies.pdf", width = 48, height = 24)

# 网络属性
attr_data <- attr_list %>% 
  map2(., names(.), \(x, y) rename(x, !!y := value)) %>% 
  reduce(\(x, y) merge(x, y, by = "name"))
write.table(attr_data, "MEN/prok.network.for_breed.genomospecies.attr.tsv", sep = "\t", row.names = F, quote = F)

plot_data <- attr_data %>% 
  filter(name != "the.number.of.keystone.nodes") %>% 
  gather(key = "group", value = "value", -name) %>% 
  left_join(select(group_x, group, class) %>% unique, by = "group") %>% 
  mutate(value = as.numeric(value))

ggbarplot(plot_data, x = "group", y = "value", fill = "group", legend = "none",
          position = position_dodge(), x.text.angle = 90, size = .4, ggtheme = theme_bw(),
          font.tickslab = c(8, "plain", "#000000"), xlab = "", ylab = "") +
  geom_hline(aes(yintercept = value),
             filter(plot_data, group == "JL-DLLW"),
             lty = "dashed", linewidth = .4) +
  facet_wrap(~factor(name), scales = "free_y") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", color = "#000000"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = .4, color = "#000000"),
        panel.background = element_rect(linewidth = .4, color = "#000000"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1/2)
ggsave("MEN/prok.network.for_breed.genomospecies.barplot.pdf", width = 15, height = 9)

# phage MEN for_breed genomospecies ------------------------

source("F:/Code/R_func/taxa.R")
source("F:/Code/R_func/calcu_MEN.R")

group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/final_vOTUs.tpm.rds")

group_x <- group %>%
  filter(group2 != "" & region == "Rectum") %>% 
  select(sample, group = breed2, class)

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

# breeds <- unique(group_x$group)
# plot_list <- list()
# attr_list <- list()
# for (i in breeds) {
#   profile_i <- select(profile_x, any_of(group_x$sample[group_x$group == i]))
#   profile_i <- profile_i[rowSums(profile_i) != 0, colSums(profile_i) != 0]
#   metadata <- data.frame(taxa = rownames(profile_i)) %>% 
#     mutate(phylum = paste0("phy_", 1:nrow(.))) %>% 
#     column_to_rownames("taxa")
#   
#   data <- calcu_MEN(profile_i, metadata, title = i) %>% suppressWarnings()
#   plot_list[[i]] <- data$gplot
#   attr_list[[i]] <- data$network_attr
# }
# save(plot_list, attr_list, file = "phage.for_breed.genomospecies.RData")

# 网络图
cowplot::plot_grid(plotlist = plot_list, nrow = 3)
ggsave("MEN/phage.network.for_breed.genomospecies.pdf", width = 48, height = 24)

# 网络属性
attr_data <- attr_list %>% 
  map2(., names(.), \(x, y) rename(x, !!y := value)) %>% 
  reduce(\(x, y) merge(x, y, by = "name"))
write.table(attr_data, "MEN/phage.network.for_breed.genomospecies.attr.tsv", sep = "\t", row.names = F, quote = F)

plot_data <- plot_data <- attr_data %>% 
  filter(name != "the.number.of.keystone.nodes") %>% 
  gather(key = "group", value = "value", -name) %>% 
  left_join(select(group_x, group, class) %>% unique, by = "group") %>% 
  mutate(value = as.numeric(value))

ggbarplot(plot_data, x = "group", y = "value", fill = "group", legend = "none",
          position = position_dodge(), x.text.angle = 90, size = .4, ggtheme = theme_bw(),
          font.tickslab = c(8, "plain", "#000000"), xlab = "", ylab = "") +
  geom_hline(aes(yintercept = value),
             filter(plot_data, group == "JL-DLLW"),
             lty = "dashed", linewidth = .4) +
  facet_wrap(~factor(name), scales = "free_y") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", color = "#000000"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = .4, color = "#000000"),
        panel.background = element_rect(linewidth = .4, color = "#000000"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1/2)
ggsave("MEN/phage.network.for_breed.genomospecies.attr.barplot.pdf", width = 15, height = 9)

# combine MEN for_breed genomospecies ------------------------

source("F:/Code/R_func/taxa.R")
source("F:/Code/R_func/calcu_MEN.R")

group <- read.delim("../profile/sample_group")
profile_x <- readRDS("../profile/genomospecies.tpm.b50.rds")
profile_y <- readRDS("../profile/final_vOTUs.tpm.rds")

group_x <- group %>%
  filter(group2 != "" & region == "Rectum") %>% 
  select(sample, group = breed2, class)

profile_x <- rbind(
  profile_x %>% 
    select(any_of(group_x$sample)) %>% 
    filter(rowSums(.) != 0),
  profile_y %>% 
    select(any_of(group_x$sample)) %>% 
    filter(rowSums(.) != 0)
)

# breeds <- unique(group_x$group)
# plot_list <- list()
# attr_list <- list()
# for (i in breeds) {
#   profile_i <- select(profile_x, any_of(group_x$sample[group_x$group == i]))
#   profile_i <- profile_i[rowSums(profile_i) != 0, colSums(profile_i) != 0]
#   metadata <- data.frame(taxa = rownames(profile_i)) %>% 
#     mutate(phylum = paste0("phy_", 1:nrow(.))) %>% 
#     column_to_rownames("taxa")
#   
#   data <- calcu_MEN(profile_i, metadata, title = i) %>% suppressWarnings()
#   plot_list[[i]] <- data$gplot
#   attr_list[[i]] <- data$network_attr
# }
# save(plot_list, attr_list, file = "combine.for_breed.genomospecies.RData")
load("MEN/combine.network.for_breed.genomospecies.RData")

# 网络图
cowplot::plot_grid(plotlist = plot_list, nrow = 3)
ggsave("MEN/combine.network.for_breed.genomospecies.pdf", width = 48, height = 24)

# 网络属性
attr_data <- attr_list %>% 
  map2(., names(.), \(x, y) rename(x, !!y := value)) %>% 
  reduce(\(x, y) merge(x, y, by = "name"))
write.table(attr_data, "MEN/combine.network.for_breed.genomospecies.attr.tsv", sep = "\t", row.names = F, quote = F)

plot_data <- gather(attr_data, key = "group", value = "value", -name) %>% 
  left_join(select(group_x, group, class) %>% unique, by = "group") %>% 
  mutate(value = ifelse(grepl("[A-Za-z]", value), 0, value)) %>% 
  mutate(value = round(as.numeric(value), 2))

ggbarplot(plot_data, x = "group", y = "value", fill = "group", legend = "none",
          position = position_dodge(), x.text.angle = 90, size = .4, ggtheme = theme_bw(),
          font.tickslab = c(8, "plain", "#000000"), xlab = "", ylab = "") +
  geom_hline(aes(yintercept = value),
             filter(plot_data, group == "JL-DLLW"),
             lty = "dashed", linewidth = .4) +
  facet_wrap(~factor(name), scales = "free_y") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", color = "#000000"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = .4, color = "#000000"),
        panel.background = element_rect(linewidth = .4, color = "#000000"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1/2)
ggsave("MEN/combine.network.for_breed.genomospecies.attr.barplot.pdf", width = 15, height = 9)

# prok MEN for_region genomospecies ------------------------

source("F:/Code/R_func/taxa.R")
source("F:/Code/R_func/calcu_MEN.R")

group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/genomospecies.tpm.b50.rds")
region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")

group_x <- group %>% 
  filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
  select(sample, group, breed = breed2, region) %>% 
  mutate(breed = factor(breed, c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")),
         region = factor(region, region_order)) %>% 
  arrange(breed, region) %>% 
  mutate(group = factor(group, unique(group)))

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

# groups <- unique(group_x$group)
# plot_list <- list()
# attr_list <- list()
# for (i in groups) {
#   profile_i <- select(profile_x, any_of(group_x$sample[group_x$group == i]))
#   profile_i <- profile_i[rowSums(profile_i) != 0, colSums(profile_i) != 0]
#   metadata <- data.frame(taxa = rownames(profile_i)) %>%
#     mutate(phylum = paste0("phy_", 1:nrow(.))) %>%
#     column_to_rownames("taxa")
# 
#   data <- calcu_MEN(profile_i, metadata, title = i) %>% suppressWarnings()
#   plot_list[[i]] <- data$gplot
#   attr_list[[i]] <- data$network_attr
# }
# save(plot_list, attr_list, file = "MEN/prok.network.for_region.genomospecies.RData")
load("MEN/prok.network.for_region.genomospecies.RData")

# 网络图
cowplot::plot_grid(plotlist = plot_list, nrow = 4)
ggsave("MEN/prok.network.for_region.genomospecies.pdf", width = 28, height = 16)

# 网络属性
attr_data <- attr_list %>% 
  map2(., names(.), \(x, y) rename(x, !!y := value)) %>% 
  reduce(\(x, y) merge(x, y, by = "name"))
write.table(attr_data, "MEN/prok.network.for_region.genomospecies.attr.tsv", sep = "\t", row.names = F, quote = F)

plot_data <- attr_data %>% 
  filter(name != "the.number.of.keystone.nodes") %>% 
  gather(key = "group", value = "value", -name) %>% 
  left_join(select(group_x, group) %>% unique, by = "group") %>% 
  mutate(value = as.numeric(value))

ggbarplot(plot_data, x = "group", y = "value", fill = "group", legend = "none",
          position = position_dodge(), x.text.angle = 90, size = .4, ggtheme = theme_bw(),
          font.tickslab = c(8, "plain", "#000000"), xlab = "", ylab = "") +
  facet_wrap(~factor(name), scales = "free_y") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", color = "#000000"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = .4, color = "#000000"),
        panel.background = element_rect(linewidth = .4, color = "#000000"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1/2)
ggsave("MEN/prok.network.for_region.genomospecies.barplot.pdf", width = 15, height = 9)

# phage MEN for_region genomospecies ------------------------

source("F:/Code/R_func/taxa.R")
source("F:/Code/R_func/calcu_MEN.R")

group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/final_vOTUs.tpm.rds")
region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")

group_x <- group %>% 
  filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
  select(sample, group, breed = breed2, region) %>% 
  mutate(breed = factor(breed, c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")),
         region = factor(region, region_order)) %>% 
  arrange(breed, region) %>% 
  mutate(group = factor(group, unique(group)))

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

# groups <- unique(group_x$group)
# plot_list <- list()
# attr_list <- list()
# for (i in groups) {
#   profile_i <- select(profile_x, any_of(group_x$sample[group_x$group == i]))
#   profile_i <- profile_i[rowSums(profile_i) != 0, colSums(profile_i) != 0]
#   metadata <- data.frame(taxa = rownames(profile_i)) %>%
#     mutate(phylum = paste0("phy_", 1:nrow(.))) %>%
#     column_to_rownames("taxa")
# 
#   data <- calcu_MEN(profile_i, metadata, title = i) %>% suppressWarnings()
#   plot_list[[i]] <- data$gplot
#   attr_list[[i]] <- data$network_attr
# }
# save(plot_list, attr_list, file = "MEN/phage.network.for_region.genomospecies.RData")
load("MEN/phage.network.for_region.genomospecies.RData")

# 网络图
cowplot::plot_grid(plotlist = plot_list, nrow = 4)
ggsave("MEN/phage.network.for_region.genomospecies.pdf", width = 28, height = 16)

# 网络属性
attr_data <- attr_list %>% 
  map2(., names(.), \(x, y) rename(x, !!y := value)) %>% 
  reduce(\(x, y) merge(x, y, by = "name"))
write.table(attr_data, "MEN/phage.network.for_region.genomospecies.attr.tsv", sep = "\t", row.names = F, quote = F)

plot_data <- plot_data <- attr_data %>% 
  filter(name != "the.number.of.keystone.nodes") %>% 
  gather(key = "group", value = "value", -name) %>% 
  left_join(select(group_x, group) %>% unique, by = "group") %>% 
  mutate(value = as.numeric(value))

ggbarplot(plot_data, x = "group", y = "value", fill = "group", legend = "none",
          position = position_dodge(), x.text.angle = 90, size = .4, ggtheme = theme_bw(),
          font.tickslab = c(8, "plain", "#000000"), xlab = "", ylab = "") +
  facet_wrap(~factor(name), scales = "free_y") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", color = "#000000"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = .4, color = "#000000"),
        panel.background = element_rect(linewidth = .4, color = "#000000"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1/2)
ggsave("MEN/phage.network.for_region.genomospecies.attr.barplot.pdf", width = 15, height = 9)

# combine MEN for_region genomospecies ------------------------

source("F:/Code/R_func/taxa.R")
source("F:/Code/R_func/calcu_MEN.R")

group <- read.delim("../profile/sample_group")
region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")

group_x <- group %>% 
  filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
  select(sample, group, breed = breed2, region) %>% 
  mutate(breed = factor(breed, c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")),
         region = factor(region, region_order)) %>% 
  arrange(breed, region) %>% 
  mutate(group = factor(group, unique(group)))

profile_x <- readRDS("../profile/genomospecies.tpm.b50.rds") %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

profile_y <- readRDS("../profile/final_vOTUs.tpm.rds") %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

samples <- intersect(colnames(profile_x), colnames(profile_y))

profile_x <- rbind(
  profile_x %>% 
    select(all_of(samples)),
  profile_y %>% 
    select(all_of(samples))
)

# groups <- unique(group_x$group)
# plot_list <- list()
# attr_list <- list()
# for (i in groups) {
#   profile_i <- select(profile_x, any_of(group_x$sample[group_x$group == i]))
#   profile_i <- profile_i[rowSums(profile_i) != 0, colSums(profile_i) != 0]
#   metadata <- data.frame(taxa = rownames(profile_i)) %>%
#     mutate(phylum = paste0("phy_", 1:nrow(.))) %>%
#     column_to_rownames("taxa")
# 
#   data <- calcu_MEN(profile_i, metadata, title = i) %>% suppressWarnings()
#   plot_list[[i]] <- data$gplot
#   attr_list[[i]] <- data$network_attr
# }
# save(plot_list, attr_list, file = "MEN/combine.network.for_region.genomospecies.RData")
load("MEN/combine.network.for_region.genomospecies.RData")

# 网络图
cowplot::plot_grid(plotlist = plot_list, nrow = 4)
ggsave("MEN/combine.network.for_region.genomospecies.pdf", width = 28, height = 24)

# 网络属性
attr_data <- attr_list %>% 
  map2(., names(.), \(x, y) rename(x, !!y := value)) %>% 
  reduce(\(x, y) merge(x, y, by = "name"))
write.table(attr_data, "MEN/combine.network.for_region.genomospecies.attr.tsv", sep = "\t", row.names = F, quote = F)

plot_data <- gather(attr_data, key = "group", value = "value", -name) %>% 
  left_join(select(group_x, group) %>% unique, by = "group") %>% 
  mutate(value = ifelse(grepl("[A-Za-z]", value), 0, value)) %>% 
  mutate(value = round(as.numeric(value), 2))

ggbarplot(plot_data, x = "group", y = "value", fill = "group", legend = "none",
          position = position_dodge(), x.text.angle = 90, size = .4, ggtheme = theme_bw(),
          font.tickslab = c(8, "plain", "#000000"), xlab = "", ylab = "") +
  facet_wrap(~factor(name), scales = "free_y") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", color = "#000000"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = .4, color = "#000000"),
        panel.background = element_rect(linewidth = .4, color = "#000000"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1/2)
ggsave("MEN/combine.network.for_region.genomospecies.attr.barplot.pdf", width = 15, height = 9)
