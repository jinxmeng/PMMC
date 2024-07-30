# Jinxin Meng, 20240301, 20240706 --------------------

setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/kegg/")
pacman::p_load(dplyr, tidyr, tibble, purrr, ggplot2, ggpubr)
library(patchwork)

kegg_db <- read.delim("/database/KEGG_v20230401/KO_level_A_B_C_D_Description", sep = "\t", quote = "")

profile <- readRDS("rawdata/KOs.tpm.rds")
group <- read.delim("../profile/sample_group")
KOs <- rownames(profile)

# KO丰度的PCoA for breed ---------------------------

source("/code/R_func/plot_PCoA.R")
source("/code/R_func/profile_process.R")
source("/code/R_func/difference_analysis.R")

group_x <- group %>%
  filter(region == "Rectum" & group2 != "") %>% 
  select(sample, group = breed2, class) %>% 
  arrange(class, group) %>% 
  mutate(group = factor(group, unique(group)))

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0) %>% 
  profile_transRA()

data <- calcu_PCoA(profile_x, group = group_x)
plot_data <- left_join(data$points, group_x, by = "sample")

colors <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#F9C851","#ACDEF3","#f881bf",
            "#94d294","#d49264","#95beda","#eb9ecd","#F9C6B3","#fff087","#c3e6dd","#b3b3b3")

p1 <- ggscatter(plot_data, x = "PCoA1", y = "PCoA2", color = "group", shape = "class",
                palette = colors, legend = "right", title = "Bray-Curtis distance-based PCoA",
                xlab = data$eig_[1], ylab = data$eig_[2], subtitle = data$label)

p2 <- ggboxplot(plot_data, x = "group", y = "PCoA1", fill = "group",
                palette = colors, rotate = T, outlier.shape = NA, color = "group",
                legend = "none") +
  stat_compare_means(ref.group = "JL-DLLW", label = "p.signif")

p3 <- ggboxplot(plot_data, x = "group", y = "PCoA2", fill = "group",
                palette = colors, x.text.angle = 90, outlier.shape = NA, color = "group",
                legend = "none") +
  stat_compare_means(ref.group = "JL-DLLW", label = "p.signif")

(p2 / p1) | (p2 / p3)

ggsave("KOs_composition/PCoA.for_breed.manual.pdf", height = 12, width = 14)

# adonis 分析

levels <- levels(group_x$group)
data <- calcu_pairwise_adonis(profile_x, group_x)

write.table(data, "KOs_composition/PCoA.for_breed.adonis.txt", sep = "\t", row.names = F, quote = F)

plot_data <- data %>% 
  mutate(x = unlist(lapply(strsplit(group_pair, "_vs_"), "[[", 1)),
         y = unlist(lapply(strsplit(group_pair, "_vs_"), "[[", 2))) %>% 
  add_plab() %>% 
  mutate(x = factor(x, levels[1:15]),
         y = factor(y, rev(levels[-1])))

ggplot(plot_data, aes(x, y)) +
  geom_tile(fill = "transparent", color = "black", width = 1, height = 1, lwd = .4) +
  geom_point(aes(size = r2adj, fill = plab), shape = 22, color = "#000000", stroke = .4) +
  scale_fill_manual(values = rev(c("#fb6a4a", "#fee0d2"))) +
  scale_size_continuous(range = c(2, 6)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1)

ggsave("KOs_composition/PCoA.for_breed.adonis.pdf", width = 6, height = 6)

# KO丰度差异 --------------------------------

source("/code/R_func/profile_process.R")
source("/code/R_func/taxa.R")
source("/code/R_func/calcu_diff.R")

group_x <- group %>%
  filter(region == "Rectum" & group2 != "") %>% 
  select(sample, group = breed2, class) %>% 
  arrange(class, group) %>% 
  mutate(group = factor(group, unique(group)))

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0) %>% 
  profile_transRA() %>% 
  mutate(lvA = kegg_db$lvA[match(rownames(.), kegg_db$lvD)]) %>% 
  group_by(lvA) %>% 
  summarise_all(sum) %>% 
  filter(!is.na(lvA))

plot_data <- profile_x %>%
  gather(key = "sample", value = "value", -lvA) %>% 
  left_join(group_x, by = "sample") %>% 
  mutate(lvAdes = kegg_db$lvAdes[match(lvA, kegg_db$lvA)]) %>% 
  filter(lvA %in% c("A09100","A09120","A09130","A09140")) %>% 
  mutate(lvAdes = factor(lvAdes, c("Metabolism","Genetic Information Processing",
                                   "Environmental Information Processing","Cellular Processes")))

ggboxplot(plot_data, x = "group", y = "value", fill = "class", width = .6, color = "class",
          palette = c("#56B4E9", "#E69F00"), outlier.shape = NA, rotate = T,
          ylab = "Relative abundance of KOs") +
  facet_wrap(vars(lvAdes), scale = "free", nrow = 1) +
  stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", hide.ns = T, color = "red")

ggsave("KOs_composition/comparison.for_breed.pdf", width = 12, height = 6)

# meta-analysis

source("/code/R_func/calcu_metafor.R")

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0) %>% 
  profile_adjacency() %>% 
  mutate(lvA = kegg_db$lvA[match(rownames(.), kegg_db$lvD)]) %>% 
  group_by(lvA) %>% 
  summarise_all(sum) %>% 
  filter(!is.na(lvA))

data <- profile_x %>% 
  column_to_rownames("lvA") %>% 
  t %>% 
  data.frame %>%
  rownames_to_column("sample") %>% 
  left_join(group_x, by = "sample") %>% 
  select(-sample, -class) %>% 
  group_by(group)
  
data2 <- data %>% 
  group_map(~cbind(
    map_dfc(.x, \(x) mean(x)) %>% 
      t %>% 
      data.frame %>% 
      rename(mean = "."),
    map_dfc(.x, \(x) sd(x)) %>% 
      t %>%
      data.frame %>% 
      rename(sd = "."), 
    map_dfc(.x, \(x) length(x)) %>% 
      t %>% 
      data.frame %>% 
      rename(n = ".")
    ) %>% 
      rownames_to_column("name")) %>% 
  map2_df(., pull(group_keys(data)), \(x, y)
       x %>% add_column(group = y)) %>% 
  ungroup() %>% 
  group_by(name)

plot_data <- data2 %>% 
  group_map(~{data.frame(filter(.x, group != "JL-DLLW") %>% 
                           select(name = group, d_mean = mean, d_sd = sd, d_n = n)) %>% 
      add_column(c_mean = filter(.x, group == "JL-DLLW") %>% pull(mean),
                 c_sd = filter(.x, group == "JL-DLLW") %>% pull(sd),
                 c_n = filter(.x, group == "JL-DLLW") %>% pull(n))} %>% 
        metafor_fit.1) %>% 
  map2_df(., pull(group_keys(data2)), \(x, y)
          x %>% add_column(index = y)) %>% 
  ungroup() %>% 
  add_plab()
  
write.table(plot_data, "KOs_composition/comparison.for_breed.meta.txt", sep = "\t", quote = F, row.names = F)

ggplot() +
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
        aspect.ratio = 1/3)

ggsave("KOs_composition/comparison.for_breed.meta.pdf", width = 8, height = 4)

# for_region by_breed 

group_x <- group %>% 
  filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
  select(sample, breed = breed2, region, class)

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0) %>% 
  profile_transRA() %>% 
  mutate(lvA = kegg_db$lvA[match(rownames(.), kegg_db$lvD)]) %>% 
  group_by(lvA) %>% 
  summarise_all(sum) %>% 
  filter(!is.na(lvA))

region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")

levels <- c("A09100","A09120","A09130","A09140")

plot_data <- profile_x %>%
  gather(key = "sample", value = "value", -lvA) %>% 
  left_join(group_x, by = "sample") %>% 
  mutate(lvAdes = kegg_db$lvAdes[match(lvA, kegg_db$lvA)],
         region = factor(region, region_order))

plot_list <- map(levels, \(x)
                 plot_data %>% 
                   filter(lvA == x) %>% 
                   ggboxplot(x = "breed", y = "value", fill = "class", width = .6, color = "class",
                             palette = c("#56B4E9", "#E69F00"), outlier.shape = NA, rotate = T,
                             ylab = "Relative abundance of KOs")  +
                   facet_grid(rows = vars(region)) +
                   stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", hide.ns = T, color = "red")
                 )

cowplot::plot_grid(plotlist = plot_list, nrow = 1)

ggsave("KOs_composition/comparison.for_region.pdf", width = 16, height = 12)

# for_region

group_x <- group %>% 
  filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
  select(sample, breed = breed2, region, class)

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0) %>% 
  profile_transRA() %>% 
  mutate(lvA = kegg_db$lvA[match(rownames(.), kegg_db$lvD)]) %>% 
  group_by(lvA) %>% 
  summarise_all(sum) %>% 
  filter(!is.na(lvA))

region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")

levels <- c("A09100","A09120","A09130","A09140")

plot_data <- profile_x %>%
  gather(key = "sample", value = "value", -lvA) %>% 
  left_join(group_x, by = "sample") %>% 
  mutate(lvAdes = kegg_db$lvAdes[match(lvA, kegg_db$lvA)],
         region = factor(region, region_order)) %>% 
  filter(lvA %in% levels)

plot_list <- map(levels, \(x) {
  comparison <- plot_data %>% 
    filter(lvA == x) %>%
    select(sample, value) %>% 
    calcu_diff(group = select(group_x, sample, group = region), group_order = region_order) %>% 
    filter(pval < 0.05) %>% 
    pull(group_pair) %>% 
    strsplit("_vs_") 
  
  plot_data %>% 
    filter(lvA == x) %>% 
    ggboxplot(x = "region", y = "value", fill = "region", width = .6, color = "region", 
              outlier.shape = NA, ylab = "Relative abundance of KOs", x.text.angle = 90, xlab = x, legend = "none",
              palette = c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#F9C851","#F9C6B3")) +
    stat_compare_means(comparisons = comparison, label = "p.signif", hide.ns = T, vjust = .8, 
                       step.increase = .05, tip.length = .01)
  }
  )

cowplot::plot_grid(plotlist = plot_list, nrow = 1)

ggsave("KOs_composition/comparison.for_region2.pdf", width = 12, height = 4)

plot_data <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0) %>% 
  colSums() %>% 
  data.frame(value = .) %>% 
  rownames_to_column("sample") %>% 
  left_join(group_x, by = "sample") %>% 
  mutate(region = factor(region, region_order))

comparison <- plot_data %>% 
  select(sample, value) %>% 
  calcu_diff(group = select(group_x, sample, group = region), group_order = region_order) %>% 
  filter(pval < 0.05) %>% 
  pull(group_pair) %>% 
  strsplit("_vs_") 

ggboxplot(plot_data, x = "region", y = "value", fill = "region", width = .6, color = "region", 
          outlier.shape = NA, ylab = "KOs TPM", x.text.angle = 90, legend = "none", xlab = "",
          palette = c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#F9C851","#F9C6B3")) +
  stat_compare_means(comparisons = comparison, label = "p.signif", hide.ns = T, vjust = .8, 
                     step.increase = .05, tip.length = .01)

ggsave("KOs_composition/comparison.for_region.tpm.pdf", width = 3, height = 4)
