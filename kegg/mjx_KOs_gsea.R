# info Jinxin Meng, 20240301, 20240705 -------------

setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/kegg/")
pacman::p_load(dplyr, tidyr, tibble, purrr, ggplot2, ggpubr)
library(ComplexHeatmap)

# GSEA 富集分析 ------------------------------------

kegg_db <- read.delim("/database/KEGG_v20230401/KO_level_A_B_C_D_Description", sep = "\t", quote = "")

files <- list.files("KOs_enrichment_gsea/", pattern = "Rectum.tsv", full.names = T)

plot_data <- data.frame(name = character(), bp = numeric(), wp = numeric())
data <- rbind()

for (i in files) {
  data_i <- read.delim(i) %>% 
    mutate(padj = p.adjust(pval, "BH"), .after = "pval") %>% 
    filter(padj < 0.05) %>% 
    mutate(enriched = ifelse(nes > 0, "bp", "wp"))
  
  if(nrow(data_i) == 0) {
    stat_i <- c(0, 0)
  } else {
    stat_i <- table(data_i$enriched) %>% 
      as.vector()
  }
  
  name <- stringr::str_split_i(i, "/", 2) %>% 
    stringr::str_split_i("_vs_", 1) %>% 
    sub("gsea_", "", x = .) %>% 
    sub("-Rectum", "", x = .)
  
  plot_data <- plot_data %>% 
    add_row(name = name, bp = stat_i[1], wp = stat_i[2])
  
  data <- rbind(data, 
                 rbind(filter(data_i, enriched == "bp") %>% 
                         select(path = Term, nes, enriched),
                       filter(data_i, enriched == "wp") %>% 
                         select(path = Term, nes, enriched)) %>% 
                   add_column(breed = name))
}

write.table(plot_data, "KOs_enrichment_gsea/diff.path.statistics.tsv", sep = "\t", quote = F, row.names = F)

plot_data <- plot_data %>% 
  arrange(desc(bp)) %>% 
  mutate(name = factor(name, name)) %>% 
  gather(key = "class", value = "value", -name) %>% 
  mutate(class = factor(class, c("bp", "wp")))

ggbarplot(plot_data, x = "name", y = "value", fill = "class", xlab = "", 
          ylab = "Number of pathway", position = position_dodge(), 
          legend = "right", legend.title = "Enriched in", label = T,
          palette = c("#56B4E9", "#E69F00"), x.text.angle = 30) +
  theme(aspect.ratio = 1/2)

ggsave("KOs_enrichment_gsea/diff.path.statistics.pdf", width = 7, height = 4)

# 富集通路 -----------------------------

data_bp <- filter(data, enriched == "bp") %>% 
  mutate(path = gsub(":.*", "", path)) %>% 
  mutate(nes = ifelse(is.infinite(nes), 8.209536, nes)) %>% 
  select(-enriched) %>% 
  spread(key = "breed", value = "nes", fill = 0) %>% 
  column_to_rownames("path")

data_bp %>% 
  rownames_to_column("ko") %>% 
  mutate(lvC = gsub("ko", "C", x = ko),
         lvB = kegg_db$lvB[match(lvC, kegg_db$lvC)],
         lvCdes = kegg_db$lvCdes[match(lvC, kegg_db$lvC)],
         lvBdes = kegg_db$lvBdes[match(lvC, kegg_db$lvC)]) %>% 
  write.table("KOs_enrichment_gsea/enriched.pathway.bp.heatmap.tsv", sep = "\t", row.names = F, quote = F)
  
data_bp <- data_bp[which(apply(data_bp, 1, \(x) sum(x != 0)) > 3), ]

pdf("KOs_enrichment_gsea/enriched.pathway.bp.heatmap2.pdf", width = 5, height = 6)
pheatmap(data_bp, scale = "none", border_color = "white", border_gp = gpar(col = "black"),
         color = circlize::colorRamp2(c(0, 8), c("white", "#56B4E9")), 
         cellheight = 6, cellwidth = 6, fontsize = 5, 
         cluster_rows = F, cluster_cols = F, 
         heatmap_legend_param = list(title = "NES value"))
dev.off()

data_wp <- filter(data, enriched == "wp") %>% 
  mutate(path = gsub(":.*", "", path)) %>% 
  mutate(nes = ifelse(is.infinite(nes), -8.209536, nes)) %>% 
  select(-enriched) %>% 
  spread(key = "breed", value = "nes", fill = 0) %>% 
  column_to_rownames("path")

data_wp %>% 
  rownames_to_column("ko") %>% 
  mutate(lvC = gsub("ko", "C", x = ko),
         lvB = kegg_db$lvB[match(lvC, kegg_db$lvC)],
         lvCdes = kegg_db$lvCdes[match(lvC, kegg_db$lvC)],
         lvBdes = kegg_db$lvBdes[match(lvC, kegg_db$lvC)]) %>% 
  write.table("KOs_enrichment_gsea/enriched.pathway.wp.heatmap.tsv", sep = "\t", row.names = F, quote = F)

data_wp <- data_wp[which(apply(data_wp, 1, \(x) sum(x != 0)) > 3), ]

pdf("KOs_enrichment_gsea/enriched.pathway.wp.heatmap2.pdf", width = 5, height = 4)
pheatmap(data_wp, scale = "none", border_color = "white", border_gp = gpar(col = "black"),
         color = circlize::colorRamp2(c(-8, 0), c("#E69F00","white")),
         cellheight = 6, cellwidth = 6, fontsize = 5, 
         cluster_rows = F, cluster_cols = F,
         heatmap_legend_param = list(title = "NES value"))
dev.off()
