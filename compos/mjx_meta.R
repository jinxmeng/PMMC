# Jinxin Meng, 20240305, 20240710 ------------------

setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/compos/")
pacman::p_load(dplyr, tidyr, tibble, purrr, ggplot2, ggpubr)
library(ComplexHeatmap)
library(ggtern)

# prok meta-analysis for_breed --------------------

source("/code/R_func/taxa.R")
source("F:/Code/R_func/calcu_metafor.R")
source("F:/Code/R_func/diversity.R")
source("F:/Code/R_func/difference_analysis.R")

group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/genomospecies.tpm.b50.rds")
taxonomy <- read.delim("../profile/genome_taxonomy.txt") %>% 
  mutate(phylum = gsub("_\\w$", "", phylum),
         family = gsub("_\\w$", "", family),
         genus = gsub("_\\w$", "", genus))
# write.table(taxonomy, "xx.tsv", sep = "\t", row.names = F, quote = F)

region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")
breed_order <- c("JL-BSWP","JL-LXP","JL-SLBP","JL-DLLW","HB-MSP","HL-MP","HL2-MP","HN-YNBP",
                 "IM-IMBP","LN-BJBP","LN-JSBP","LN-LDBP","LN-THNP","LN2-LDBP","SD-LWP","SX-LLBP")

group_x <- group %>%
  filter(group2 != "" & region == "Rectum") %>% 
  select(sample, group = breed2)

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0) %>% 
  profile_transRA()

# phylum|family|genus|genomospecies 

taxa <- "phylum"
ref <- "JL-DLLW"
breeds <- setdiff(breed_order, ref)

data <- taxa_trans(profile_x, taxonomy, to = taxa, out_all = T, transRA = T)

meta <- t(data) %>% 
  data.frame(check.names = F) %>% 
  mutate(group = group_x$group[match(rownames(.), group_x$sample)]) %>% 
  group_by(group) %>% 
  group_modify(~map2_df(.x, colnames(.x), \(x, y) 
                        data.frame(name = y, mean = mean(x), sd = sd(x), n = length(x)))) %>% 
  ungroup

name <- unique(meta$name)
meta_out <- rbind()
for (i in name) {
  meta_ref <- meta %>% 
    filter(name == i, group %in% ref) %>% 
    select(mean, sd, n) %>% 
    unlist(use.names = F)
  meta_data <- meta %>% 
    filter(name == i, group %in% breeds) %>% 
    rename(d_mean = mean, d_sd = sd, d_n = n) %>% 
    add_column(c_mean = meta_ref[1], c_sd = meta_ref[2], c_n = meta_ref[3])
  
  smd_meta <- escalc(measure = "SMD", data = meta_data, append = T, 
                     m1i = d_mean, m2i = c_mean, 
                     sd1i = d_sd, sd2i = c_sd, 
                     n1i = d_n, n2i = c_n) %>% 
    na.omit %>% 
    data.frame
  
  non_na <- smd_meta %>% 
    filter(!is.na(yi)) %>% 
    nrow()
  
  if(non_na != 0){
    smd_rma <- rma(yi, vi, method = "REML", data = smd_meta, control = list(stepadj = 0.5, maxiter = 10000))
    # merge each data
    smd_meta <- smd_rma$data %>% 
      add_column(measure = "SMD", # 效应值的计算方法
                 model = "RM",  # 累积效应值计算模型
                 method_tau2 = "REML", # 随机效应模型估计案例内方差（Tau^2）的计算方法
                 val_tau2 = as.numeric(smd_rma$tau2), # 案例间方差的值
                 I2 = paste0(round(smd_rma$I2, digits = 2), "%"), # 案例间差异大小占总差异的比例
                 Q = smd_rma$QE, # 异质性检验
                 Q_pval = smd_rma$QEp, # 异质性检验P值 越显著异质性越大
                 estimate = as.numeric(smd_rma$beta),
                 ci_lb = smd_rma$ci.lb,
                 ci_ub = smd_rma$ci.ub,
                 pval = smd_rma$pval)
    meta_out <- rbind(meta_out, smd_meta)
  }
}

meta_out <- meta_out %>% 
  mutate(padj = p.adjust(pval, method = "BH")) %>% 
  add_plab(by = "padj")
write.table(meta_out, paste0("meta_analysis/meta.for_breed.", taxa, "_smd.tsv"), sep = "\t", quote = F, row.names = F)

plot_data <- meta_out %>%
  select(group, name, yi, estimate, ci_lb, ci_ub, padj, plab) %>%
  filter(!is.na(yi))
  # filter(plab != "")

ggplot() +
  geom_vline(xintercept = 0, color = "#000000", lwd = .4, lty = 2) +
  geom_point(aes(x = estimate, y = name),
             plot_data %>% select(name, estimate) %>% unique(),
             size = 2.5, shape = 15, inherit.aes = F) +
  geom_errorbar(aes(y = name, xmin = ci_lb, xmax = ci_ub),
                plot_data %>% select(name, ci_lb, ci_ub) %>% unique(),
                width = .2, lwd = .4, inherit.aes = F) +
  geom_point(aes(x = yi, y = name, color = group),
             plot_data,
             size = 2, inherit.aes = F, show.legend = F) +
  geom_text(aes(x = min(yi), y = name, label =  plab),
            plot_data, color = "red") +
  labs(y = "", x = "Standradized Mean Difference (Random Effect Model)") +
  scale_color_viridis_d(begin = .4) +
  theme_classic2() +
  theme(axis.text = element_text(color = "black", size = 8),
        axis.title = element_text(size = 8),
        axis.line = element_line(linewidth = .4, color = "#000000"),
        axis.ticks = element_line(linewidth = .4, color = "#000000"))
ggsave(paste0("meta_analysis/meta.for_breed.", taxa, "_smd.pdf"), width = 6,
       height = length(name)/5, limitsize = FALSE)

# phage meta-analysis for_breed --------------------

source("/code/R_func/taxa.R")
source("F:/Code/R_func/calcu_metafor.R")
source("F:/Code/R_func/difference_analysis.R")

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

# family|genomospecies 

taxa = "name"
ref <- "JL-DLLW"
breeds <- setdiff(breed_order, ref)

# data <- taxa_trans(profile_x, taxonomy, to = taxa, out_all = T, transRA = T)
data <- profile_x

meta <- t(data) %>% 
  data.frame(check.names = F) %>% 
  mutate(group = group_x$group[match(rownames(.), group_x$sample)]) %>% 
  group_by(group) %>% 
  group_modify(~map2_df(.x, colnames(.x), \(x, y) 
                        data.frame(name = y, mean = mean(x), sd = sd(x), n = length(x)))) %>% 
  ungroup

name <- unique(meta$name)
meta_out <- list()
pb <- txtProgressBar(style = 3)
x = 0
for (i in name) {
  x = x + 1
  meta_ref <- meta %>% 
    filter(name == i, group %in% ref) %>% 
    select(mean, sd, n) %>% 
    unlist(use.names = F)
  meta_data <- meta %>% 
    filter(name == i, group %in% breeds) %>% 
    rename(d_mean = mean, d_sd = sd, d_n = n) %>% 
    add_column(c_mean = meta_ref[1], c_sd = meta_ref[2], c_n = meta_ref[3])
  
  smd_meta <- escalc(measure = "SMD", data = meta_data, append = T, 
                     m1i = d_mean, m2i = c_mean, 
                     sd1i = d_sd, sd2i = c_sd, 
                     n1i = d_n, n2i = c_n) %>% 
    na.omit %>% 
    data.frame
  
  non_na <- smd_meta %>% 
    filter(!is.na(yi)) %>% 
    nrow()
  
  if(non_na != 0){
    smd_rma <- rma(yi, vi, method = "REML", data = smd_meta, control = list(stepadj = 0.5, maxiter = 10000))
    # merge each data
    smd_meta <- smd_rma$data %>% 
      add_column(measure = "SMD", # 效应值的计算方法
                 model = "RM",  # 累积效应值计算模型
                 method_tau2 = "REML", # 随机效应模型估计案例内方差（Tau^2）的计算方法
                 val_tau2 = as.numeric(smd_rma$tau2), # 案例间方差的值
                 I2 = paste0(round(smd_rma$I2, digits = 2), "%"), # 案例间差异大小占总差异的比例
                 Q = smd_rma$QE, # 异质性检验
                 Q_pval = smd_rma$QEp, # 异质性检验P值 越显著异质性越大
                 estimate = as.numeric(smd_rma$beta),
                 ci_lb = smd_rma$ci.lb,
                 ci_ub = smd_rma$ci.ub,
                 pval = smd_rma$pval)
    meta_out[[i]] <- smd_meta
    setTxtProgressBar(pb, x/length(name))
  }
}
close(pb)

meta_out <- map_dfr(meta_out, \(x) x)
meta_out <- meta_out %>% 
  mutate(padj = p.adjust(pval, method = "BH")) %>% 
  add_plab(by = "padj")
write.table(meta_out, paste0("meta_analysis/phage.meta.for_breed.", taxa, "_smd.tsv"), sep = "\t", quote = F, row.names = F)

# prok meta stat -----------------------

taxa_level <- c("phylum", "family", "genus", "genomospecies")

data <- map(taxa_level, \(x) 
    a = read.delim(paste0("meta_analysis/meta.for_breed.", x, "_smd.tsv")) %>% 
      select(name, estimate, padj) %>% 
      filter(padj < 0.05) %>% 
      distinct() %>% 
      mutate(enriched = ifelse(estimate > 0, "NBP", "CP")) %>% 
      pull(enriched) %>% 
      table()
    )

names(data) <- taxa_level

taxa_total <- taxonomy %>% 
  select(all_of(taxa_level)) %>% 
  map_df(\(x) length(unique(x))) %>% 
  t %>% 
  data.frame(total = .) %>% 
  rownames_to_column("taxa")

map_df(data, \(x) x) %>% 
  add_column(taxa = taxa_level) %>% 
  left_join(taxa_total, by = "taxa") %>% 
  rowwise() %>% 
  mutate(none = total - CP - NBP) %>% 
  select(-total) %>% 
  gather(key = "class", value = "value", -taxa) %>% 
  mutate(class = factor(class, c("NBP", "CP", "none")),
         taxa = factor(taxa, taxa_level)) %>% 
  ggbarplot("class", "value", fill = "class", palette = c("#56B4E9", "#E69F00", "grey"), 
            position = position_dodge(), label = T, xlab = "", ylab = "Number of features") +
  facet_wrap(vars(taxa), scale = "free_y", nrow = 1) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 10))
ggsave("meta_analysis/meta.stat.pdf", height = 4, width = 8)

# phage meta stat -----------------------
taxa_level <- c("family", "name")

data <- map(taxa_level, \(x) 
            read.delim(paste0("meta_analysis/phage.meta.for_breed.", x, "_smd.tsv")) %>% 
              select(name, estimate, padj) %>% 
              filter(padj < 0.05) %>% 
              distinct() %>% 
              mutate(enriched = ifelse(estimate > 0, "NBP", "CP")) %>% 
              pull(enriched) %>% 
              table()
)

names(data) <- taxa_level

taxa_total <- taxonomy %>% 
  select(all_of(taxa_level)) %>% 
  map_df(\(x) length(unique(x))) %>% 
  t %>% 
  data.frame(total = .) %>% 
  rownames_to_column("taxa")

map_df(data, \(x) x) %>% 
  add_column(taxa = taxa_level) %>% 
  left_join(taxa_total, by = "taxa") %>% 
  rowwise() %>% 
  mutate(none = total - CP - NBP) %>% 
  select(-total) %>% 
  gather(key = "class", value = "value", -taxa) %>% 
  mutate(class = factor(class, c("NBP", "CP", "none")),
         taxa = factor(taxa, taxa_level)) %>% 
  ggbarplot("class", "value", fill = "class", palette = c("#56B4E9", "#E69F00", "grey"), 
            position = position_dodge(), label = T, xlab = "", ylab = "Number of features") +
  facet_wrap(vars(taxa), scale = "free_y", nrow = 1) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 10))
ggsave("meta_analysis/phage.meta.stat.pdf", height = 4, width = 5)

# prok meta genus volcano ----------------------

data <- read.delim("meta_analysis/meta.for_breed.genus_smd.tsv")
taxonomy <- read.delim("../profile/genome_taxonomy.txt") %>% 
  mutate(phylum = gsub("_\\w$", "", phylum),
         family = gsub("_\\w$", "", family),
         genus = gsub("_\\w$", "", genus))

isolate <- taxonomy %>% 
  mutate(source = ifelse(grepl("PN", name), "SGBs", "isolates")) %>% 
  group_by(genus) %>% 
  group_modify(~.x %>% 
              mutate(label = ifelse("isolates" %in% source, "yes", "no"))) %>% 
  filter(label == "yes") %>% 
  pull(genus) %>% 
  unique()
  
plot_data <- data %>% 
  select(estimate, name, padj) %>% 
  distinct() %>% 
  mutate(enriched = ifelse(padj > 0.05, "none", ifelse(estimate > 0, "NBP", "CP"))) %>% 
  mutate(estimate = ifelse(estimate < -3.5, -3.5, estimate),
         padj = ifelse(-log10(padj) > 50, 1e-50, padj),
         label = ifelse(name %in% isolate & padj < 0.05, name, ""))

colors <- structure(c("#56B4E9", "#E69F00", "grey"), names = c("NBP", "CP", "none"))
labels <- table(plot_data$enriched) %>% 
  data.frame() %>% 
  mutate(label = paste0(Var1, ", ", Freq)) %>% 
  pull(label, name = Var1)

ggplot(plot_data, aes(estimate, -log10(padj), color = enriched)) +
  geom_point(size = 2, alpha = .8) +
  geom_text(aes(label = label), color = "black", size = 1) +
  scale_color_manual(values = colors, labels = labels) +
  scale_y_continuous(expand = c(.01, .01)) +
  labs(x = "Estimate of Meta-analysis", y = "-Log10-transformed P.adj", color = "Enriched in") + 
  geom_hline(yintercept = -log10(0.05), lty = 2, lwd = .4, color = "black") +
  theme_pubr() +
  theme(aspect.ratio = 0.9)
ggsave("meta_analysis/meta.for_breed.genus.volcano.pdf", height = 6, width = 7)

taxonomy_x <- taxonomy %>% 
  filter(genus %in% c("g__Turicibacter", "g__Akkermansia", "g__Lachnoclostridium"))

group_x <- group %>%
  filter(group2 != "" & region == "Rectum") %>% 
  select(sample, group = breed2)

profile_x <- readRDS("../profile/genomospecies.tpm.rds") %>% 
  filter(rownames(.) %in% taxonomy_x$genomospecies) %>% 
  select(all_of(group_x$sample)) %>% 
  profile_smp2grp(group_x)

diff_label <- readRDS("../profile/genomospecies.tpm.rds") %>% 
  filter(rownames(.) %in% taxonomy_x$genomospecies) %>% 
  select(all_of(group_x$sample)) %>% 
  calcu_diff_profile(., group_x, center_group = "JL-DLLW") %>%
  mutate(group = stringr::str_split_i(group_pair, "_vs_", 1)) %>% 
  add_plab(format = 3) %>% 
  select(group, plab, name) %>% 
  spread(key = "group", value = "plab") %>% 
  add_column(`JL-DLLW` = "") %>% 
  column_to_rownames("name")

diff_label <- diff_label[rownames(profile_x), colnames(profile_x)]

annotation_row <- taxonomy_x %>% 
  select(genomospecies, genus) %>% 
  mutate(genomospecies = factor(genomospecies, rownames(diff_label))) %>% 
  arrange(genomospecies) %>% 
  column_to_rownames("genomospecies")
split_row <- factor(annotation_row$genus)

pdf("meta_analysis/meta.for_breed.genus.species.ab2.pdf", width = 8, height = 8)
pheatmap(t(profile_x), scale = "column",
         color = colorRampPalette(c("#74add1", "#ffffff", "#f46d43"))(100),
         show_rownames = T, show_colnames = T, 
         fontsize_row = 8, fontface_col = 8, 
         show_row_dend = F, show_column_dend = F,
         cellheight = 10, cellwidth = 10,
         column_gap = unit(0, "mm"),
         column_split = split_row,
         border_gp = gpar(col = "black"),
         border_color = "white",
         row_title_rot = 0,
         cluster_column_slices = F,
         display_numbers = as.matrix(t(diff_label)))
dev.off()

# core bac --------------------------

source("/code/R_func/taxa.R")
source("/code/R_func/calcu_diff.R")

group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/genomospecies.tpm.rds")
taxonomy <- read.delim("../profile/genome_taxonomy.txt") %>% 
  mutate(phylum = gsub("_\\w$", "", phylum),
         family = gsub("_\\w$", "", family))

region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")
breed_order <- c("JL-BSWP","JL-LXP","JL-SLBP","JL-DLLW","HB-MSP","HL-MP","HL2-MP","HN-YNBP",
                 "IM-IMBP","LN-BJBP","LN-JSBP","LN-LDBP","LN-THNP","LN2-LDBP","SD-LWP","SX-LLBP")

group_x <- group %>%
  filter(region == "Rectum" & group2 != "") %>% 
  select(sample, group = breed2, class) %>% 
  arrange(class, group) %>% 
  mutate(group = factor(group, unique(group)))

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

x <- profile_prevalence(profile, by_group = F) %>% left_join(taxonomy, by = "name")

profile_x["G66",] %>%
  t %>%
  data.frame() %>%
  rownames_to_column("sample") %>% 
  left_join(group_x, "sample") %>% 
  ggboxplot("group", "G66", fill = "group") +
  stat_compare_means(ref.group = "JL-DLLW", label = "p.signif")



