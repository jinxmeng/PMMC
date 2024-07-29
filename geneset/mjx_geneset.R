# Jinxin Meng, 20240503, 20240702 -------------------------------

pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr)
setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/geneset/")

# 基因集合聚类比较 --------------------------------------

# 在服务器上跑的
LiM_2024 <- read.delim("LiM_2024.gene.list", header = F)
HuJ_2024 <- read.delim("HuJ_2024.gene.list", header = F)
ChenC_2021 <- read.delim("ChenC_2021.gene.list", header = F)
XiaoL_2016 <- read.delim("XiaoL_2016.gene.list", header = F)

genes <- list(LiM_2024 = LiM_2024$V1,
              XiaoL_2016 = XiaoL_2016$V1,
              ChenC_2021 = ChenC_2021$V1,
              HuJ_2024 = HuJ_2024$V1)

library(ggVennDiagram)

p <- ggVennDiagram(genes, label = "both", label_color = "black", label_alpha = 0,edge_lty = "dashed", edge_size = 1) +
  scale_fill_gradient(low="#fef5f5",high = "#f0ad80", name = "Number of gene clusters")
ggsave(p, "geneset.contrast.venn.pdf", width = 5, height = 5)

# 稀释曲线分析 样本水平 -----------------------------

read.delim("rawdata/gene.subsample.tsv") %>% 
  mutate(obs = obs/1000000, sd = sd/1000000) %>% 
  ggline(x = "step", y = "obs", point.size = NA,
       xlab = "Number of samples", ylab = "Number of non-redundant genes (×10^6)")

ggsave("rare/gene.subsample.pdf", width = 4, height = 4)

# 物种注释饼图 -----------------------------------------

source("/Code/R_func/plot_pie.R")

plot_data <- read.delim("geneset.tax.txt")

plot_data %>% 
  select(name, n) %>% 
  plot_pie(fill = "auto", add_n = T, font_size = 4)

ggsave("geneset.tax.pdf", width = 4, height = 4)

# 功能注释柱状图 ------------------------------------------

colors <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#F9C851","#F9C6B3")

plot_data <- read.delim("geneset.function.txt") %>% 
  arrange(desc(n))

ggbarplot(plot_data, x = "name", y = "n_million", fill = "name", legend = "none", label = plot_data$prec,
          ylab = "Number of gene clusters (×106)", xlab = "Functional gene categroy") +
  scale_fill_manual(values = colors <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#F9C851","#F9C6B3")) +
  theme(aspect.ratio = 3/4)

ggsave("geneset.function.pdf", width = 4, height = 3)

# 组装和基因预测统计 ------------------------------------

group <- read.delim("../profile/sample_group")

plot_data <- read.delim("geneset.metadata.txt") %>% 
  left_join(., group, by = "sample")

p1 <- ggscatter(plot_data, x = "contigs_avg_len", y = "contigs_N50_len", fill = "breed2", shape = 21,
          size = 2.5, xlab = "Average length of contigs (bp)", ylab = "N50 length of contigs (bp)") +
  theme_gray() +
  theme(aspect.ratio = 3/2,
        axis.text = element_text(color = "black"))

p2 <- ggscatter(plot_data, x = "contigs", y = "gene", fill = "breed2", shape = 21,
                size = 2.5, xlab = "Number of contigs", ylab = "Number of predicted genes") +
  theme_gray() +
  theme(aspect.ratio = 3/2,
        axis.text = element_text(color = "black"))

p3 <- ggscatter(plot_data, x = "contigs_avg_len", y = "gene_avg_len", fill = "breed2", shape = 21,
                size = 2.5, xlab = "Average length of contigs (bp)", ylab = "Average length of predicted genes (bp)") +
  theme_gray() +
  theme(aspect.ratio = 3/2,
        axis.text = element_text(color = "black"))

p4 <- ggscatter(plot_data, x = "contigs_avg_len", y = "contigs_N50_len", fill = "region", shape = 21,
                size = 2.5, xlab = "Average length of contigs (bp)", ylab = "N50 length of contigs (bp)") +
  theme_gray() +
  theme(aspect.ratio = 3/2,
        axis.text = element_text(color = "black"))

p5 <- ggscatter(plot_data, x = "contigs", y = "gene", fill = "region", shape = 21,
                size = 2.5, xlab = "Number of contigs", ylab = "Number of predicted genes") +
  theme_gray() +
  theme(aspect.ratio = 3/2,
        axis.text = element_text(color = "black"))

p6 <- ggscatter(plot_data, x = "contigs_avg_len", y = "gene_avg_len", fill = "region", shape = 21,
                size = 2.5, xlab = "Average length of contigs (bp)", ylab = "Average length of predicted genes (bp)") +
  theme_gray() +
  theme(aspect.ratio = 3/2,
        axis.text = element_text(color = "black"))

cowplot::plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2)
ggsave("geneset.contigs.gene.scatter.pdf", width = 14, height = 8)

# 基因在各组的丰富度 ------------------------------------

source("/code/R_func/calcu_metafor.R")
source("/code/R_func/difference_analysis.R")

group <- read.delim("../profile/sample_group")

# 不同品种
group_x <- group %>%
  filter(region == "Rectum" & group2 != "") %>% 
  select(sample, group = breed2, class) %>% 
  arrange(class, group) %>% 
  mutate(group = factor(group, unique(group)))

metadata <- read.delim("geneset.metadata.txt") %>% 
  select(sample, gene_clusters) %>% 
  filter(sample %in% group_x$sample)

data <- left_join(metadata, group_x, by = "sample") %>% 
  mutate(n = gene_clusters/1e6)

var <- data %>% 
  filter(grepl("PN5R", sample)) %>% 
  pull(n) %>% 
  mean()

ggboxplot(data, x = "group", y = "n", xlab = "", ylab = "Richness (×10^6)",
          fill = "class", color = "class", legend = "none", x.text.angle = 30, 
          palette = c("#56B4E9", "#E69F00"), outlier.shape = NA) +
  stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", hide.ns = T, color = "red") +
  geom_hline(yintercept = var, lty = "dashed", linewidth = .4) +
  theme(aspect.ratio = 1/2)

ggsave("richness/comparison.for_breed.pdf", width = 6, height = 4)

# meta分析

vec <- data %>% 
  filter(group == "JL-DLLW") %>% 
  pull(n) 

plot_data <- data %>% 
  select(-class, -sample) %>% 
  filter(group != "JL-DLLW") %>% 
  group_by(group) %>% 
  summarise(d_mean = mean(n), d_sd = sd(n), d_n = n()) %>% 
  ungroup %>% 
  add_column(c_mean = mean(vec),
             c_sd = sd(vec),
             c_n = length(vec)) %>% 
  rename(name = group) %>% 
  metafor_fit.1() %>% 
  add_plab() %>% 
  add_column(index = "richness")

write.table(plot_data, "richness/comparison.for_breed.meta.txt", sep = "\t", row.names = F, quote = F)

ggplot(plot_data) +
  geom_vline(xintercept = 0, color = "#000000", lwd = .4, lty = 2) +
  geom_point(aes(x = estimate, y = index), 
             plot_data %>% select(index, estimate) %>% unique(), 
             size = 2.5, shape = 15, inherit.aes = F) + 
  geom_errorbar(aes(y = index, xmin = ci_lb, xmax = ci_ub),
                plot_data %>% select(index, ci_lb, ci_ub) %>% unique(),
                width = .2, lwd = .4, inherit.aes = F) +
  geom_point(aes(x = yi, y = index, color = name),
             plot_data, 
             size = 2, inherit.aes = F, show.legend = F) +
  geom_text(aes(x = min(yi), y = index, label =  plab), 
            plot_data) +
  labs(y = "", x = "Standradized Mean Difference (Random Effect Model)") +
  scale_color_viridis_d(begin = .4) +
  theme_classic2() +
  theme(axis.text = element_text(color = "black", size = 8),
        axis.title = element_text(size = 8),
        axis.line = element_line(linewidth = .4, color = "#000000"),
        axis.ticks = element_line(linewidth = .4, color = "#000000"),
        aspect.ratio = 1/6)

ggsave("richness/comparison.for_breed.meta.pdf", height = 2, width = 6)

# 不同肠段

group_x <- group %>% 
  filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
  select(sample, breed = breed2, region, class)

metadata <- read.delim("geneset.metadata.txt") %>% 
  select(sample, gene_clusters) %>% 
  filter(sample %in% group_x$sample)

region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")

data <- left_join(metadata, group_x, by = "sample") %>% 
  mutate(n = gene_clusters/1e6,
         region = factor(region, region_order))

plot_data <- data %>% 
  group_by(breed, region) %>% 
  summarise(mean = mean(n), sd = sd(n), n = n()) %>% 
  mutate(se = sd / sqrt(n))

ggerrorplot(plot_data, x = "region", y = "mean", color = "breed", 
                  position = position_dodge(width = 0), xlab = "Location", 
                  ylab = "Richness (×10^6)", x.text.angle = 90) +
  geom_ribbon(aes(x = region, ymin = mean - se, ymax = mean + se, fill = breed, group = breed),
              plot_data, show.legend = F, inherit.aes = F, alpha = .3) +
  stat_compare_means(aes(x = region, y = n, group = breed, color = breed), 
                     data, inherit.aes = F,
                     label = "p.signif") +
  scale_color_manual(values = c("#82ccdc","#b4b3d8","#c2b75f","#F9C851")) +
  scale_fill_manual(values = c("#82ccdc","#b4b3d8","#c2b75f","#F9C851")) +
  theme(aspect.ratio = 1)

ggsave("richness/comparison.for_region.pdf", width = 6, height = 6)

ggboxplot(data, x = "breed", y = "n", ylab = "Richness (×10^6)", 
          fill = "class", color = "class", legend = "none", x.text.angle = 30, 
          palette = c("#56B4E9", "#E69F00"), outlier.shape = NA) +
  facet_wrap(vars(region), scales = "free_x", nrow = 1) +
  stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", hide.ns = T, color = "red")

ggsave("richness/comparison.for_region2.pdf", width = 9, height = 5)

# meta分析

data_c <- data %>% 
  filter(breed == "JL-DLLW") %>%
  group_by(breed, region) %>% 
  summarise(c_mean = mean(n), c_sd = sd(n), c_n = n()) %>% 
  ungroup %>% 
  rbind(., ., .)

plot_data <- data %>% 
  select(-class, -sample) %>% 
  filter(breed != "JL-DLLW") %>% 
  group_by(breed, region) %>% 
  summarise(d_mean = mean(n), d_sd = sd(n), d_n = n()) %>% 
  ungroup %>% 
  cbind(select(data_c, starts_with("c_"))) %>% 
  arrange(region) %>% 
  group_by(region) %>% 
  group_modify(~metafor_fit.1(rename(.x, name = breed))) %>% 
  ungroup() %>% 
  rename(breed = name) %>% 
  add_plab()
  
write.table(plot_data, "richness/comparison.for_region.meta.txt", sep = "\t", row.names = F, quote = F)

ggplot(plot_data) +
  geom_vline(xintercept = 0, color = "#000000", lwd = .4, lty = 2) +
  geom_point(aes(x = estimate, y = region), 
             plot_data %>% select(region, estimate) %>% unique(), 
             size = 2.5, shape = 15, inherit.aes = F) + 
  geom_errorbar(aes(y = region, xmin = ci_lb, xmax = ci_ub),
                plot_data %>% select(region, ci_lb, ci_ub) %>% unique(),
                width = .2, lwd = .4, inherit.aes = F) +
  geom_point(aes(x = yi, y = region, color = breed),
             plot_data, 
             size = 2, inherit.aes = F, show.legend = F) +
  geom_text(aes(x = min(yi), y = region, label =  plab), 
            plot_data) +
  labs(y = "Region", x = "Standradized Mean Difference (Random Effect Model)") +
  scale_color_viridis_d(begin = .4) +
  theme_classic2() +
  theme(axis.text = element_text(color = "black", size = 8),
        axis.title = element_text(size = 8),
        axis.line = element_line(linewidth = .4, color = "#000000"),
        axis.ticks = element_line(linewidth = .4, color = "#000000"),
        aspect.ratio = 1/2)

ggsave("richness/comparison.for_region.meta.pdf", height = 4, width = 6)
