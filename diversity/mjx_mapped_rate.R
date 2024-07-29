# Jinxin Meng, 20240626, 20240709 -----------------------

setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/diversity/")
pacman::p_load(dplyr, tidyr, tibble, purrr, ggplot2, ggpubr, ggpmisc)

# mapped rate ------------------------

group <- read.delim("../profile/sample_group")
group_x <- filter(group, class == "white", region == "Rectum")

# prok

data <- read.delim("rawdata/prok.stat.mapped.rate.tsv", header = F, 
                   col.names = c("sample", "mapped_reads", "total_reads", "mapped_rate")) %>%
  filter(sample %in% group_x$sample) %>% 
  mutate(mapped_rate = mapped_rate * 100) %>% 
  left_join(group_x, by = "sample")

ggboxplot(data, x = "group", y = "mapped_rate", fill = "group", legend = "none",
          palette = c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f"), 
          ylab = "Mapped rate (%)", xlab = "", x.text.angle = 90) +
  stat_compare_means(ref.group = "JL-DLLW-Rectum", label = "p.signif", color = "red") +
  theme(aspect.ratio = 1)

ggsave("mapped_rate/prok.mapped_rate.pdf", height = 5, width = 6)

# phage

data <- read.delim("rawdata/phage.stat.mapped.rate.tsv", header = F, 
                   col.names = c("sample", "mapped_reads", "total_reads", "mapped_rate")) %>%
  filter(sample %in% group_x$sample) %>% 
  mutate(mapped_rate = mapped_rate * 100) %>% 
  left_join(group_x, by = "sample")

ggboxplot(data, x = "group", y = "mapped_rate", fill = "group", legend = "none",
          palette = c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f"), 
          ylab = "Mapped rate (%)", xlab = "", x.text.angle = 90) +
  stat_compare_means(ref.group = "JL-DLLW-Rectum", label = "p.signif", color = "red") +
  theme(aspect.ratio = 1)

ggsave("mapped_rate/phage.mapped_rate.bar.pdf", height = 5, width = 6)

