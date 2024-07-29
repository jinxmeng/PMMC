# Jinxin Meng, 20240305, 20240713 ------------------

setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/compos/")
pacman::p_load(dplyr, tidyr, tibble, purrr, ggplot2, ggpubr)
library(ComplexHeatmap)
library(ggtern)
library(scatterpie)

# pork diff test for breed --------------------

source("/code/R_func/taxa.R")
source("/code/R_func/calcu_diff.R")

group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/genomospecies.tpm.b50.rds")
taxonomy <- read.delim("../profile/genome_taxonomy.txt") %>% 
  mutate(phylum = gsub("_\\w$", "", phylum),
         family = gsub("_\\w$", "", family),
         genus = gsub("_\\w$", "", genus))

region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")
breed_order <- c("JL-BSWP","JL-LXP","JL-SLBP","JL-DLLW","HB-MSP","HL-MP","HL2-MP","HN-YNBP",
                 "IM-IMBP","LN-BJBP","LN-JSBP","LN-LDBP","LN-THNP","LN2-LDBP","SD-LWP","SX-LLBP")

group_x <- group %>%
  filter(group2 != "" & region == "Rectum") %>% 
  select(sample, group = breed2)

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

# phylum level ----------------------

data <- taxa_trans(profile_x, taxonomy, to = "phylum", out_all = T, transRA = T) %>% 
  profile_filter(group_x, by_group = T, n_group = 2, min_n = 3)

diff <- calcu_diff_profile(data, group_x, center_group = "JL-DLLW", plab = T, padj = T)

write.table(diff, "diff_test/prok.diff_test.phylum.tsv", sep = "\t", quote = F, row.names = F)

plot_data <- profile_smp2grp(data, group_x)

plot_plab <- mutate(diff, group_pair = unlist(lapply(strsplit(group_pair, "_vs_"), "[[", 1))) %>% 
  select(-pval, -padj, -method) %>% 
  spread(key = "group_pair", value = "plab", fill = "") %>% 
  column_to_rownames(var = "name") %>% 
  add_column(`JL-DLLW` = "") %>% 
  data.frame(check.names = F)
plot_plab <- plot_plab[rownames(plot_data), colnames(plot_data)]
plot_plab %>% 
  rownames_to_column("name") %>% 
  write.table("diff_test/prok.diff_test.phylum.plab.tsv", sep = "\t", quote = F, row.names = F)

vec_ref <- plot_data["JL-DLLW"] %>% 
  unlist(use.names = F)
sign_data <- apply(plot_data, 2, \(x) 
                   ifelse((x - vec_ref) > 0, "+", 
                          ifelse((x - vec_ref) < 0, "-", ""))) %>% 
  data.frame(check.names = F)

plot_bar <- apply(plot_plab, 2, \(x) ifelse(x!="", "*", "")) %>% 
  data.frame(check.names = F) %>% 
  map2_df(sign_data, ., \(x, y) paste0(x, y)) %>% 
  data.frame(row.names = rownames(sign_data), check.names = F) %>% 
  t() %>%
  data.frame() %>% 
  map2_df(., colnames(.), \(x, y) data.frame(name = y, up = sum(x=="+*"), down = sum(x=="-*")) %>% 
            mutate(none = 15 - up - down))
write.table(plot_bar, "diff_test/prok.diff_test.phylum.stat.txt", sep = "\t", row.names = F, quote = F) 

plot_bar %>%
  gather(key = "enriched", value = "value", -name) %>%
  mutate(enriched = factor(enriched, c("up", "down", "none") %>% rev)) %>% 
  ggbarplot("name", "value", fill = "enriched", x.text.angle = 60,
            palette = c("#56B4E9", "#E69F00", "grey") %>% rev)
ggsave("diff_test/prok.diff_test.phylum.stat.pdf", width = 5, height = 4)

anno_bar <- columnAnnotation(bar = anno_barplot(column_to_rownames(plot_bar, var = "name"), gp = gpar(fill = c("#56B4E9", "#E69F00", "grey"))))
anno_bar_lgd <- Legend(at = c("NBP", "CP", "none"), legend_gp = gpar(fill = c("#56B4E9", "#E69F00", "grey")), title = "stat")

pdf("diff_test/prok.diff_test.phylum.heatmap.pdf", width = 6, height = 6)
pheatmap(t(plot_data), scale = "column",
         border_color = "black", 
         color = colorRampPalette(c("#74add1", "#ffffff", "#f46d43"))(100),
         border_gp = gpar(col = "black"),
         cellwidth = 15, cellheight = 15, 
         show_row_dend = F, show_column_dend = F,
         display_numbers = t(plot_plab),
         top_annotation = anno_bar,
         fontsize = 10, fontsize_col = 8, fontsize_row = 8,
         heatmap_legend_param = list(title = "Scale relative abundance (%)",
                                     title_gp = grid::gpar(fontface = "plain"),
                                     title_position = "leftcenter-rot"))
draw(anno_bar_lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))
dev.off()

# family level ----------------------

data <- taxa_trans(profile_x, taxonomy, to = "family", out_all = T, transRA = T) %>% 
  profile_filter(group_x, by_group = T, n_group = 2, min_n = 3)

diff <- calcu_diff_profile(data, group_x, center_group = "JL-DLLW", plab = T, padj = T)

write.table(diff, "diff_test/prok.diff_test.family.tsv", sep = "\t", quote = F, row.names = F)

plot_data <- profile_smp2grp(data, group_x)

plot_plab <- mutate(diff, group_pair = unlist(lapply(strsplit(group_pair, "_vs_"), "[[", 1))) %>% 
  select(-pval, -padj, -method) %>% 
  spread(key = "group_pair", value = "plab", fill = "") %>% 
  column_to_rownames(var = "name") %>% 
  add_column(`JL-DLLW` = "") %>% 
  data.frame(check.names = F)
plot_plab <- plot_plab[rownames(plot_data), colnames(plot_data)]
plot_plab %>% 
  rownames_to_column("name") %>% 
  write.table("diff_test/prok.diff_test.family.plab.tsv", sep = "\t", quote = F, row.names = F)

vec_ref <- plot_data["JL-DLLW"] %>% 
  unlist(use.names = F)
sign_data <- apply(plot_data, 2, \(x) 
                   ifelse((x - vec_ref) > 0, "+", 
                          ifelse((x - vec_ref) < 0, "-", ""))) %>% 
  data.frame(check.names = F)

plot_bar <- apply(plot_plab, 2, \(x) ifelse(x!="", "*", "")) %>% 
  data.frame(check.names = F) %>% 
  map2_df(sign_data, ., \(x, y) paste0(x, y)) %>% 
  data.frame(row.names = rownames(sign_data), check.names = F) %>% 
  t() %>%
  data.frame(check.names = F) %>% 
  map2_df(., colnames(.), \(x, y) data.frame(name = y, up = sum(x=="+*"), down = sum(x=="-*")) %>% 
        mutate(none = 15 - up - down))
write.table(plot_bar, "diff_test/prok.diff_test.family.stat.txt", sep = "\t", row.names = F, quote = F) 

plot_bar %>%
  gather(key = "enriched", value = "value", -name) %>%
  mutate(enriched = factor(enriched, c("up", "down", "none") %>% rev)) %>% 
  ggbarplot("name", "value", fill = "enriched", x.text.angle = 90,
            palette = c("#56B4E9", "#E69F00", "grey") %>% rev)
ggsave("diff_test/prok.diff_test.family.stat.pdf", width = 20, height = 7)

anno_bar <- columnAnnotation(bar = anno_barplot(column_to_rownames(plot_bar, var = "name"), gp = gpar(fill = c("#56B4E9", "#E69F00", "grey"))))
anno_bar_lgd <- Legend(at = c("NBP", "CP", "none"), legend_gp = gpar(fill = c("#56B4E9", "#E69F00", "grey")), title = "stat")

pdf("diff_test/prok.diff_test.family.heatmap.pdf", width = 18, height = 6)
pheatmap(t(plot_data), scale = "column",
         border_color = "black", 
         color = colorRampPalette(c("#74add1", "#ffffff", "#f46d43"))(100),
         border_gp = gpar(col = "black"),
         cellwidth = 8, cellheight = 8, 
         show_row_dend = F, show_column_dend = F,
         display_numbers = t(plot_plab),
         top_annotation = anno_bar,
         fontsize = 6, fontsize_col = 6, fontsize_row = 6,
         heatmap_legend_param = list(title = "Scale relative abundance (%)",
                                     title_gp = grid::gpar(fontface = "plain"),
                                     title_position = "leftcenter-rot"))
draw(anno_bar_lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))
dev.off()

plot_bar %>% 
  group_by(up, down) %>% 
  summarise(n = n()) %>% 
  data.frame() %>% 
  write.table("diff_test/prok.diff_test.family.enriched_pattern.txt", sep = "\t", row.names = F, quote = F)

taxonomy_x <- taxonomy %>% 
  select(phylum, family) %>%
  mutate(phylum = forcats::fct_lump_n(phylum, 9, other_level = "p__Other")) %>% 
  unique()

taxa_order <- levels(taxonomy_x$phylum)

plot_data <- plot_bar %>% 
  left_join(taxonomy_x, by = c("name" = "family")) %>% 
  group_by(up, down, none) %>% 
  group_modify(~ cbind(data.frame(n = nrow(.x)), 
              table(.x$phylum) %>% 
              data.frame() %>% 
              column_to_rownames("Var1") %>%
              t %>% 
              data.frame(row.names = NULL))
              ) %>% 
  ungroup

colors <- c("#bebada","#ffffb3","#8dd3c7","#FF9D9A","#fee08b",
            "#a7c9e0","#de77ae","#fabfd2","#b8e186","#d9d9d9")

ggplot() + 
  geom_scatterpie(aes(up, down, r = scales::rescale(n, to = c(.2, .5))), plot_data, cols = taxa_order) + 
  scale_fill_manual(values = colors) +
  labs(x = "up", y = "down") +
  geom_scatterpie_legend(scales::rescale(plot_data$n, to = c(.2, .5)), n = 4, x = 14, y = 14) +
  coord_fixed() +
  theme_pubr(legend = "right")
ggsave("diff_test/prok.diff_test.family.enriched_pattern.pdf", width = 8, height = 7)

plot_data <- plot_bar %>% 
  group_by(up, down) %>% 
  summarise(n = n()) %>% 
  ungroup %>%
  rowwise() %>% 
  mutate(x = paste0("up_", up, "&down_", down))

ggbarplot(plot_data, x = "x", y = "n", fill = "x", legend = "none", x.text.angle = 90, 
          xlab = "pattern", ylab = "Number of features", label =  T)

ggtern(plot_data, aes(up, down, none)) +
  geom_point(aes(fill = phylum), shape = 21, color = "#000000", stroke = .4) +
  scale_fill_manual(values = colors) +
  scale_size_continuous(range = c(1, 6)) +
  theme_rgbw()


# genus level ----------------------

data <- taxa_trans(profile_x, taxonomy, to = "genus", out_all = T, transRA = T) %>% 
  profile_filter(group_x, by_group = T, n_group = 2, min_n = 3)

diff <- calcu_diff_profile(data, group_x, center_group = "JL-DLLW", plab = T, padj = T)

write.table(diff, "diff_test/prok.diff_test.genus.tsv", sep = "\t", quote = F, row.names = F)

plot_data <- profile_smp2grp(data, group_x)

plot_plab <- mutate(diff, group_pair = unlist(lapply(strsplit(group_pair, "_vs_"), "[[", 1))) %>% 
  select(-pval, -padj, -method) %>% 
  spread(key = "group_pair", value = "plab", fill = "") %>% 
  column_to_rownames(var = "name") %>% 
  add_column(`JL-DLLW` = "") %>% 
  data.frame(check.names = F)
plot_plab <- plot_plab[rownames(plot_data), colnames(plot_data)]
plot_plab %>% 
  rownames_to_column("name") %>% 
  write.table("diff_test/prok.diff_test.genus.plab.tsv", sep = "\t", quote = F, row.names = F)

vec_ref <- plot_data["JL-DLLW"] %>% 
  unlist(use.names = F)
sign_data <- apply(plot_data, 2, \(x) 
                   ifelse((x - vec_ref) > 0, "+", 
                          ifelse((x - vec_ref) < 0, "-", ""))) %>% 
  data.frame(check.names = F)

stat_data <- apply(plot_plab, 2, \(x) ifelse(x!="", "*", "")) %>% 
  data.frame(check.names = F) %>% 
  map2_df(sign_data, ., \(x, y) paste0(x, y)) %>% 
  data.frame(row.names = rownames(sign_data), check.names = F) %>% 
  t() %>%
  data.frame(check.names = F) %>% 
  map2_df(., colnames(.), \(x, y) data.frame(name = y, up = sum(x=="+*"), down = sum(x=="-*")) %>% 
            mutate(none = 15 - up - down))
write.table(stat_data, "diff_test/prok.diff_test.genus.stat.txt", sep = "\t", row.names = F, quote = F) 

stat_data %>% 
  group_by(up, down) %>% 
  summarise(n = n()) %>% 
  data.frame() %>% 
  write.table("diff_test/prok.diff_test.genus.enriched_pattern.txt", sep = "\t", row.names = F, quote = F)

taxonomy_x <- taxonomy %>% 
  select(phylum, genus) %>%
  mutate(phylum = forcats::fct_lump_n(phylum, 9, other_level = "p__Other")) %>% 
  unique()

taxa_order <- levels(taxonomy_x$phylum)

plot_data <- stat_data %>% 
  left_join(taxonomy_x, by = c("name" = "genus")) %>% 
  group_by(up, down, none) %>% 
  group_modify(~ cbind(data.frame(n = nrow(.x)), 
                       table(.x$phylum) %>% 
                         data.frame() %>% 
                         column_to_rownames("Var1") %>%
                         t %>% 
                         data.frame(row.names = NULL))
  ) %>% 
  ungroup

colors <- c("#bebada","#ffffb3","#8dd3c7","#FF9D9A","#fee08b",
            "#a7c9e0","#de77ae","#fabfd2","#b8e186","#d9d9d9")

ggplot() + 
  geom_scatterpie(aes(up, down, r = scales::rescale(n, to = c(.2, .5))), plot_data, cols = taxa_order) + 
  scale_fill_manual(values = colors) +
  labs(x = "up", y = "down") +
  geom_scatterpie_legend(scales::rescale(plot_data$n, to = c(.2, .5)), n = 4, x = 14, y = 14) +
  coord_fixed() +
  theme_pubr(legend = "right")
ggsave("diff_test/prok.diff_test.genus.enriched_pattern.pdf", width = 8, height = 7)

# genomospecies level ----------------------

data <- taxa_trans(profile_x, taxonomy, to = "genomospecies", out_all = T, transRA = T) %>% 
  profile_filter(group_x, by_group = T, n_group = 2, min_n = 3)

diff <- calcu_diff_profile(data, group_x, center_group = "JL-DLLW", plab = T, padj = T)

write.table(diff, "diff_test/prok.diff_test.genomospecies.tsv", sep = "\t", quote = F, row.names = F)

plot_data <- profile_smp2grp(data, group_x)

plot_plab <- mutate(diff, group_pair = unlist(lapply(strsplit(group_pair, "_vs_"), "[[", 1))) %>% 
  select(-pval, -padj, -method) %>% 
  spread(key = "group_pair", value = "plab", fill = "") %>% 
  column_to_rownames(var = "name") %>% 
  add_column(`JL-DLLW` = "") %>% 
  data.frame(check.names = F)
plot_plab <- plot_plab[rownames(plot_data), colnames(plot_data)]
plot_plab %>% 
  rownames_to_column("name") %>% 
  write.table("diff_test/prok.diff_test.genomospecies.plab.tsv", sep = "\t", quote = F, row.names = F)

vec_ref <- plot_data["JL-DLLW"] %>% 
  unlist(use.names = F)
sign_data <- apply(plot_data, 2, \(x) 
                   ifelse((x - vec_ref) > 0, "+", 
                          ifelse((x - vec_ref) < 0, "-", ""))) %>% 
  data.frame(check.names = F)

stat_data <- apply(plot_plab, 2, \(x) ifelse(x!="", "*", "")) %>% 
  data.frame(check.names = F) %>% 
  map2_df(sign_data, ., \(x, y) paste0(x, y)) %>% 
  data.frame(row.names = rownames(sign_data), check.names = F) %>% 
  t() %>%
  data.frame(check.names = F) %>% 
  map2_df(., colnames(.), \(x, y) data.frame(name = y, up = sum(x=="+*"), down = sum(x=="-*")) %>% 
            mutate(none = 15 - up - down))
write.table(stat_data, "diff_test/prok.diff_test.genomospecies.stat.txt", sep = "\t", row.names = F, quote = F) 

stat_data %>% 
  group_by(up, down) %>% 
  summarise(n = n()) %>% 
  data.frame() %>% 
  write.table("diff_test/prok.diff_test.genomospecies.enriched_pattern.txt", sep = "\t", row.names = F, quote = F)

taxonomy_x <- taxonomy %>% 
  select(phylum, genomospecies) %>%
  mutate(phylum = forcats::fct_lump_n(phylum, 9, other_level = "p__Other")) %>% 
  unique()

taxa_order <- levels(taxonomy_x$phylum)

plot_data <- stat_data %>% 
  left_join(taxonomy_x, by = c("name" = "genomospecies")) %>% 
  group_by(up, down, none) %>% 
  group_modify(~ cbind(data.frame(n = nrow(.x)), 
                       table(.x$phylum) %>% 
                         data.frame() %>% 
                         column_to_rownames("Var1") %>%
                         t %>% 
                         data.frame(row.names = NULL))
  ) %>% 
  ungroup

colors <- c("#bebada","#ffffb3","#8dd3c7","#FF9D9A","#fee08b",
            "#a7c9e0","#de77ae","#fabfd2","#b8e186","#d9d9d9")

ggplot() + 
  geom_scatterpie(aes(up, down, r = scales::rescale(n, to = c(.2, .5))), plot_data, cols = taxa_order) + 
  scale_fill_manual(values = colors) +
  labs(x = "up", y = "down") +
  geom_scatterpie_legend(scales::rescale(plot_data$n, to = c(.2, .5)), n = 4, x = 14, y = 14) +
  coord_fixed() +
  theme_pubr(legend = "right")
ggsave("diff_test/prok.diff_test.genomospecies.enriched_pattern.pdf", width = 8, height = 7)

# phage diff test for breed --------------------

source("/code/R_func/taxa.R")
source("/code/R_func/calcu_diff.R")

group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/final_vOTUs.tpm.rds")
taxonomy <- read.delim("../phage_genome/data.phage_genome_metadata.txt") %>% 
  select(name, family)

region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")
breed_order <- c("JL-BSWP","JL-LXP","JL-SLBP","JL-DLLW","HB-MSP","HL-MP","HL2-MP","HN-YNBP",
                 "IM-IMBP","LN-BJBP","LN-JSBP","LN-LDBP","LN-THNP","LN2-LDBP","SD-LWP","SX-LLBP")

group_x <- group %>%
  filter(group2 != "" & region == "Rectum") %>% 
  select(sample, group = breed2)

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

# family level ----------------------

data <- taxa_trans(profile_x, taxonomy, to = "family", out_all = T, transRA = T) %>% 
  profile_filter(group_x, by_group = T, n_group = 2, min_n = 3)

diff <- calcu_diff_profile(data, group_x, center_group = "JL-DLLW", plab = T, padj = T)

write.table(diff, "diff_test/phage.diff_test.family.tsv", sep = "\t", quote = F, row.names = F)

plot_data <- profile_smp2grp(data, group_x)

plot_plab <- mutate(diff, group_pair = unlist(lapply(strsplit(group_pair, "_vs_"), "[[", 1))) %>% 
  select(-pval, -padj, -method) %>% 
  spread(key = "group_pair", value = "plab", fill = "") %>% 
  column_to_rownames(var = "name") %>% 
  add_column(`JL-DLLW` = "") %>% 
  data.frame(check.names = F)
plot_plab <- plot_plab[rownames(plot_data), colnames(plot_data)]
plot_plab %>% 
  rownames_to_column("name") %>% 
  write.table("diff_test/phage.diff_test.family.plab.tsv", sep = "\t", quote = F, row.names = F)

vec_ref <- plot_data["JL-DLLW"] %>% 
  unlist(use.names = F)
sign_data <- apply(plot_data, 2, \(x) 
                   ifelse((x - vec_ref) > 0, "+", 
                          ifelse((x - vec_ref) < 0, "-", ""))) %>% 
  data.frame(check.names = F)

plot_bar <- apply(plot_plab, 2, \(x) ifelse(x!="", "*", "")) %>% 
  data.frame(check.names = F) %>% 
  map2_df(sign_data, ., \(x, y) paste0(x, y)) %>% 
  data.frame(row.names = rownames(sign_data), check.names = F) %>% 
  t() %>%
  data.frame(check.names = F) %>% 
  map2_df(., colnames(.), \(x, y) data.frame(name = y, up = sum(x=="+*"), down = sum(x=="-*")) %>% 
            mutate(none = 15 - up - down))
write.table(plot_bar, "diff_test/phage.diff_test.family.stat.txt", sep = "\t", row.names = F, quote = F) 

plot_bar %>%
  gather(key = "enriched", value = "value", -name) %>%
  mutate(enriched = factor(enriched, c("up", "down", "none") %>% rev)) %>% 
  ggbarplot("name", "value", fill = "enriched", x.text.angle = 90,
            palette = c("#56B4E9", "#E69F00", "grey") %>% rev)
ggsave("diff_test/phage.diff_test.family.stat.pdf", width = 8, height = 7)

anno_bar <- columnAnnotation(bar = anno_barplot(column_to_rownames(plot_bar, var = "name"), gp = gpar(fill = c("#56B4E9", "#E69F00", "grey"))))
anno_bar_lgd <- Legend(at = c("NBP", "CP", "none"), legend_gp = gpar(fill = c("#56B4E9", "#E69F00", "grey")), title = "stat")

pdf("diff_test/phage.diff_test.family.heatmap.pdf", width = 7, height = 6)
pheatmap(t(plot_data), scale = "column",
         border_color = "black", 
         color = colorRampPalette(c("#74add1", "#ffffff", "#f46d43"))(100),
         border_gp = gpar(col = "black"),
         cellwidth = 15, cellheight = 15, 
         show_row_dend = F, show_column_dend = F,
         display_numbers = t(plot_plab),
         top_annotation = anno_bar,
         fontsize = 10, fontsize_col = 8, fontsize_row = 8,
         heatmap_legend_param = list(title = "Scale relative abundance (%)",
                                     title_gp = grid::gpar(fontface = "plain"),
                                     title_position = "leftcenter-rot"))
draw(anno_bar_lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))
dev.off()

# genomospecies level ----------------------

data <- profile_filter(profile_x, group_x, by_group = T, n_group = 2, min_n = 3)

diff <- calcu_diff_profile(data, group_x, center_group = "JL-DLLW", plab = T, padj = T)

write.table(diff, "diff_test/phage.diff_test.genomospecies.tsv", sep = "\t", quote = F, row.names = F)

plot_data <- profile_smp2grp(data, group_x)

plot_plab <- mutate(diff, group_pair = unlist(lapply(strsplit(group_pair, "_vs_"), "[[", 1))) %>% 
  select(-pval, -padj, -method) %>% 
  spread(key = "group_pair", value = "plab", fill = "") %>% 
  column_to_rownames(var = "name") %>% 
  add_column(`JL-DLLW` = "") %>% 
  data.frame(check.names = F)
plot_plab <- plot_plab[rownames(plot_data), colnames(plot_data)]
plot_plab %>% 
  rownames_to_column("name") %>% 
  write.table("diff_test/phage.diff_test.genomospecies.plab.tsv", sep = "\t", quote = F, row.names = F)

vec_ref <- plot_data["JL-DLLW"] %>% 
  unlist(use.names = F)
sign_data <- apply(plot_data, 2, \(x) 
                   ifelse((x - vec_ref) > 0, "+", 
                          ifelse((x - vec_ref) < 0, "-", ""))) %>% 
  data.frame(check.names = F)

plot_bar <- apply(plot_plab, 2, \(x) ifelse(x!="", "*", "")) %>% 
  data.frame(check.names = F) %>% 
  map2_df(sign_data, ., \(x, y) paste0(x, y)) %>% 
  data.frame(row.names = rownames(sign_data), check.names = F) %>% 
  t() %>%
  data.frame(check.names = F) %>% 
  map2_df(., colnames(.), \(x, y) data.frame(name = y, up = sum(x=="+*"), down = sum(x=="-*")) %>% 
            mutate(none = 15 - up - down))
write.table(plot_bar, "diff_test/phage.diff_test.genomospecies.stat.txt", sep = "\t", row.names = F, quote = F) 
