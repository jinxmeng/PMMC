# Jinxin Meng, 20240324, 20240628 -----------------------------

pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr)
setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/phage_genome/")
source("/code/R_func/plot_pie.R")

# 病毒质量 --------------------------------------------------

data <- read.delim("data.phage_genome_metadata.txt")

data$checkv_quality %>% 
  table(name = .) %>% 
  data.frame %>% 
  rename(n = Freq) %>% 
  plot_pie(color = "black", font_size = 4, fill = c("#82ccdc","#f0ad80","#b4b3d8"), start = 90)

ggsave("phage.vOTUs.qs.pdf", width = 4, height = 4)

# 聚类比较 PVD -----------------------------------------------

# 根据聚类结果模拟一个韦恩图
# genomospecies
# /share/data1/mjx/proj/04.black_pig_metagenome_20230529/20.vir_cluster/compare_to_PVD/total_vOTUs.cls
# PVD-specific vOTUs: 45919 
# Bpig-specific vOTUs: 8625
# shared vOTUs: 2252

library(ggvenn)

x = list("LiM_2024" = 1:10877, "PVD" = 8626:56796)

ggvenn(x, fill_color = c("#F9C6B3","#ACDEF3"), auto_scale = T, stroke_color = "white", stroke_size = .8, )

ggsave("phage.cluster_contrast.venn.pdf", width = 5, height = 5)

# 病毒物种 ------------------------------------

data <- read.delim("data.phage_genome_metadata.txt")

data$family %>% 
  table(name = .) %>% 
  data.frame() %>% 
  ggbarplot(x = "name", y = "Freq", fill = "name", label = T, size = .5,
            legend = "none", xlab = "", ylab = "Number of vOTUs",
            x.text.angle = 90, sort.by.groups = F, sort.val = "desc")

ggsave("phage.number_vOTUs_at_family.bar.pdf", width = 8, height = 5)

# 病毒宿主 ------------------------------------

data <- read.delim("data.phage_genome_metadata.txt") %>% 
  select(family, host) %>% 
  mutate(host = gsub("_\\w$", "", host),
         family = forcats::fct_lump_n(family, 9, ties.method = "first"),
         host = forcats::fct_lump_n(host, 10, ties.method = "first")) %>% 
  group_by(family, host) %>% 
  summarise(value = n()) %>% 
  ungroup()

taxa_order <- data %>% 
  select(family, value) %>% 
  group_by(family) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  arrange(desc(value)) %>% 
  pull(family)

host_order <- data %>% 
  select(host, value) %>% 
  group_by(host) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  arrange(desc(value)) %>% 
  pull(host) %>% 
  as.character()

colors <- c("#A0CBE8","#FFBE7D","#8CD17D","#F1CE63","#86BCB6","#ffff99",
            "#FF9D9A","#FABFD2","#D4A6C8","#D7B5A6","#b3b3b3")

data %>% 
  mutate(family = factor(family, taxa_order),
         host = factor(host, host_order)) %>% 
  ggbarplot(x = "family", y = "value", fill = "host", legend = "right", x.text.angle= 30,
            palette = colors, xlab = "", ylab = "Number of virus")

ggsave("phage.taxa_host.pdf", width = 7, height = 4.5)


# 病毒类型 ------------------------------------

data$type %>% 
  table(name = .) %>% 
  data.frame() %>%
  plot_pie(dat_colnames = c("n" = "Freq"), fill = c("#82ccdc","#f0ad80"),
           color = "black", add_n = T, font_size = 4)

ggsave("phage.vOTUs.type.pdf", height = 4, width = 4)

# 病毒AMG ------------------------------------

kegg_db <- read.delim("/database/KEGG_v20230401/KO_level_A_B_C_D_Description", sep = "\t", quote = "")
metadata <- read.delim("data.phage_genome_metadata.txt")

data <- read.delim("rawdata/VIBRANT_AMG_individuals_virus.cls.tsv") %>% 
  filter(scaffold %in% metadata$name) 
write.table(data, "phage.AMGs.table.tsv", sep = "\t", quote = F, row.names = F)

data <- data %>% 
  select(name = scaffold, ko = AMG.KO) %>% 
  mutate(family = metadata$family[match(name, metadata$name)],
         lv2des = kegg_db$lvBdes[match(ko, kegg_db$lvD)],
         name = metadata$vOTU_rename[match(name, metadata$name)])
  
data_n <- data$name %>% 
  table(name = .) %>% 
  data.frame() %>% 
  arrange(desc(Freq)) %>% 
  filter(Freq >= 8)

plot_data <- data %>%
  filter(name %in% data_n$name) %>% 
  mutate(name = paste0(name, " (", family, ")")) %>% 
  group_by(name, lv2des) %>% 
  summarise(n = n()) %>% 
  spread(key = "lv2des", value = "n", fill = 0) %>% 
  ungroup() %>% 
  column_to_rownames("name")

library(ComplexHeatmap)
pdf("phage.AMGs.table.phage_gt8.pdf", width = 8, height = 10)
pheatmap(plot_data, 
         breaks = c(-2, 1, 4, 7, 10, 13, 16),
         color = c("#ffffff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6"),
         border_gp = gpar(col = "black"),
         border_color = "#ffffff",
         cluster_rows = F, cluster_cols = F, 
         heatmap_legend_param = list(title = "Number of AMGs"),
         cellwidth = 10, cellheight = 10)
dev.off()

data$lv2des %>% 
  table %>% 
  data.frame %>% 
  rename(name = ".") %>% 
  arrange(desc(Freq)) %>% 
  mutate(name = factor(name, name)) %>% 
  ggbarplot("name", "Freq", x.text.angle = 60, legend = "none",
            sort.val = "desc", sort.by.groups = T, fill = "name",
            xlab = "", ylab = "Number of AMGs")

ggsave("phage.AMGs.function.pdf", width = 8, height = 6)

data <- data %>% 
  mutate(family = forcats::fct_lump_n(family, 6),
         lv2des = forcats::fct_lump_n(lv2des, 7))

data_n <- data %>% 
  select(name, family) %>% 
  unique() %>% 
  group_by(family) %>% 
  summarise(value = n()) %>% 
  ungroup() %>% 
  arrange(desc(value))

plot_data <- data %>% 
  left_join(data_n, by = "family") %>% 
  mutate(family = paste0(family, " (n=", value, ")"))

function_order <- plot_data %>% 
  select(lv2des) %>%
  table %>% 
  sort(decreasing = T) %>% 
  names()

plot_data %>% 
  group_by(family, lv2des) %>% 
  summarise(value = n()) %>% 
  ungroup() %>% 
  mutate(family = factor(family, data_n %>% 
                           mutate(family = paste0(family = paste0(family, " (n=", value, ")"))) %>% 
                           pull(family)),
         lv2des = factor(lv2des, function_order)) %>% 
  ggbarplot("family", "value", fill = "lv2des", x.text.angle = 60, legend = "right",
            palette = c("#A0CBE8","#FFBE7D","#8CD17D","#F1CE63","#86BCB6","#FF9D9A","#FABFD2","#D4A6C8"),
            xlab = "", ylab = "Number of KOs")

ggsave("phage.AMGs.taxa.function.pdf", width = 8, height = 7)


# 病毒重命名 -------------------------

metadata <- read.delim("data.phage_genome_metadata.txt") %>% 
  mutate(vOTU_rename = paste0("vOTU_", sprintf("%05d", 1:nrow(.))))
write.table(metadata, "data.phage_genome_metadata.txt", sep = "\t", row.names = F, quote = F)

# 进化树 ----------------------------
# 
# library(ggtree)
# library(ggtreeExtra)
# library(tidytree)
# library(ggnewscale)
# library(treeio)
# 
# tr <- read.tree("../profile/final_vOTUs.all.bionj.desc.newick")
# metadata <- read.delim("data.phage_genome_metadata.txt") %>% 
#   mutate(name = gsub("\\|", "_", name))
# 
# tr_bar_family <- select(metadata, name, family) %>% 
#   mutate(family = fct_lump_n(family, 9, other_level = "f__Other"))
# tr_family_colors <- tr_bar_family$family %>% 
#   table %>% 
#   sort(decreasing = T) %>%
#   names() %>% 
#   discard(~ . == "f__Other") %>% 
#   append("f__Other") %>% 
#   structure(c("#58BBD0","#EB9256","#9B99CB","#AD9F2A","#87CCBA","#ffed6f",
#               "#7fc97f","#e78ac3","#80b1d3","#d9d9d9"), names = .)
# 
# tr_bar_type <- select(metadata, name, type)
# tr_type_colors <- structure(c("#f1948e","#4dbbd6"), names = c("lytic", "lysogenic"))
# 
# tr_bar_expand <- select(metadata, genome, expand)
# tr_expand_colors <- structure(c("#66c2a5","#ffffff"), names = c("Yes", "No"))
# 
# p <- ggtree(tr, layout = "fan", lwd = .4, branch.length = "none") +
#   geom_rootedge(rootedge = .06) +
#   geom_fruit(aes(y = name, fill = family),
#              tr_bar_family, 
#              geom_tile, 
#              pwidth = 3, offset = .1, inherit.aes = F) +
#   scale_fill_manual(values = tr_family_colors) +
#   new_scale_fill() +
#   geom_fruit(aes(y = name, fill = type),
#              tr_bar_type,
#              geom_tile, 
#              pwidth = 3, offset = .1, inherit.aes = F) +
#   scale_fill_manual(values = tr_type_colors) +
#   new_scale_fill()
#   geom_fruit(aes(y = genome, fill = expand),
#              tr_bar_expand,
#              geom_tile,
#              pwidth = .09, offset = .15, inherit.aes = F) +
#   scale_fill_manual(values = tr_expand_colors) +
#   new_scale_fill() +
#   geom_fruit(aes(y = genome, fill = novo),
#              tr_bar_novo,
#              geom_tile,
#              pwidth = .09, offset = .15, inherit.aes = F) +
#   scale_fill_manual(values = tr_novo_colors) +
#   new_scale_fill() +
#   geom_fruit(aes(y = genome, x = n),
#              tr_bar_count, 
#              geom_bar, 
#              stat = "identity", orientation = "y", pwidth = .4, offset = .15, 
#              inherit.aes = F, grid.params = list(color = "#56b4e9", size = .6))
# 
# open_tree(p, angle = 180)
# ggsave("phage.ggtree.v1.pdf", width = 20, height = 20)  

