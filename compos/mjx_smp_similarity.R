# Jinxin Meng, 20240305, 20240628 -------------------------------------

setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/compos/")
pacman::p_load(dplyr, tidyr, tibble, purrr, ggplot2, ggpubr)

# prok 层次聚类 ----------------------------------------

source("/code/R_func/taxa.R")
group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/genomospecies.tpm.b50.rds")
taxonomy <- read.delim("../profile/genome_taxonomy.txt")

# hclust
library(ggtree)

plot_list <- list()
for (i in colnames(taxonomy)[-c(1,2)]) {
  plot_list[[i]] <- taxa_trans(profile, taxonomy, group, to = i, out_all = T, smp2grp = T, transRA = T) %>% 
    t %>% 
    dist %>% 
    hclust %>% 
    treeio::as.phylo(.) %>% 
    ggtree(layout = "fan", linetype = "dashed", color = "#487AA1") +
    geom_tippoint(aes(fill = label), color = "#000000", size = 2, shape = 21, show.legend = F) +
    geom_tiplab(align = T, hjust = 0, offset = 2, size = 3, linetype = NA) + 
    labs(caption = i) + 
    theme(aspect.ratio = 1, 
          plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "cm"))
}

cowplot::plot_grid(plotlist = plot_list, nrow = 2)
ggsave("prok.smps_similarity.hclust_tr.pdf", width = 15, height = 10)

# prok 相关性 -------------------------------------------------

library(linkET)
plot_list <- list()
for (i in colnames(taxonomy)[-c(1,2)]) {
  plot_list[[i]] <- taxa_trans(profile, taxonomy, group, to = i, out_all = T, smp2grp = T, transRA = T) %>% 
    correlate() %>% 
    qcorrplot(type = "lower", diag = FALSE) + 
    geom_point(aes(size = r, fill = r), shape = 21) + 
    scale_fill_viridis_c() + 
    scale_size_continuous(range = c(1,4)) + 
    labs(caption = i) + 
    theme(aspect.ratio = 1)
}
cowplot::plot_grid(plotlist = plot_list, nrow = 2)
ggsave("prok.smps_similarity.corr.pdf", width = 30, height = 15)

# prok 聚类热图 ---------------------------------------------

library(ComplexHeatmap)
plot_list <- list()
for (i in colnames(taxonomy)[-c(1,2)]) {
  plot_list[[i]] <- taxa_trans(profile, taxonomy, group, to = i, top_n = 100, smp2grp = T, transRA = T) %>% 
    pheatmap(scale = "row", show_rownames = F, cellwidth = 6, cellheight = 3, fontsize_col = 5,
             treeheight_col = 20, treeheight_row = 20, main = i, border_color = NA,
             heatmap_legend_param = list(title = ""),) %>% 
    draw() %>% 
    grid.grabExpr()
}

pdf("prok.smps_similarity.cluster_heatmap_top100.pdf", width = 14, height = 12)
gridExtra::grid.arrange(grobs = plot_list, nrow = 2, ncol = 3)
dev.off()

# phage 层次聚类 ----------------------------------------

source("/code/R_func/taxa.R")
group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/final_vOTUs.tpm.rds")
taxonomy <- read.delim("../phage_genome/data.phage_genome_metadata.txt") %>% 
  select(name, family)

library(ggtree)

p1 <- taxa_trans(profile, taxonomy, group, to = "family", out_all = T, smp2grp = T, transRA = T) %>% 
  t %>% 
  dist %>% 
  hclust %>% 
  treeio::as.phylo(.) %>% 
  ggtree(layout = "fan", linetype = "dashed", color = "#487AA1") +
  geom_tippoint(aes(fill = label), color = "#000000", size = 2, shape = 21, show.legend = F) +
  geom_tiplab(align = T, hjust = 0, offset = 2, size = 3, linetype = NA) + 
  labs(caption = "family") + 
  theme(aspect.ratio = 1, 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "cm"))

p2 <- profile_smp2grp(profile, group) %>% 
  profile_transRA() %>% 
  t %>% 
  dist %>% 
  hclust %>% 
  treeio::as.phylo(.) %>% 
  ggtree(layout = "fan", linetype = "dashed", color = "#487AA1") +
  geom_tippoint(aes(fill = label), color = "#000000", size = 2, shape = 21, show.legend = F) +
  geom_tiplab(align = T, hjust = 0, offset = 2, size = 3, linetype = NA) + 
  labs(caption = "vOTUs") + 
  theme(aspect.ratio = 1, 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "cm"))

cowplot::plot_grid(p1, p2, nrow = 1)

ggsave("phage.smps_similarity.hclust_tr.pdf", width = 10, height = 5)

# phage 相关性 -------------------------------------------------

library(linkET)

p1 <- taxa_trans(profile, taxonomy, group, to = "family", out_all = T, smp2grp = T, transRA = T) %>% 
  correlate() %>% 
  qcorrplot(type = "lower", diag = FALSE) + 
  geom_point(aes(size = r, fill = r), shape = 21) + 
  scale_fill_viridis_c() + 
  scale_size_continuous(range = c(1,4)) + 
  labs(caption = "family") + 
  theme(aspect.ratio = 1)

ggsave(p1, filename = "phage.smps_similarity.corr.pdf", width = 10, height = 10)

# phage 聚类热图 ---------------------------------------------

library(ComplexHeatmap)

pdf("phage.smps_similarity.cluster_heatmap_top100.pdf", width = 8, height = 4)
taxa_trans(profile, taxonomy, group, to = "family", top_n = 100, smp2grp = T, transRA = T) %>% 
    pheatmap(scale = "row", show_rownames = F, cellwidth = 6, cellheight = 3, fontsize_col = 5,
             treeheight_col = 20, treeheight_row = 20, main = "family", border_color = NA,
             heatmap_legend_param = list(title = ""))
dev.off()
