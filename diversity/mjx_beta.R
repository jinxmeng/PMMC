# Jinxin Meng, 20240626, 20240709 --------------------

setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/diversity/")
pacman::p_load(dplyr, tidyr, tibble, purrr, ggplot2, ggpubr, ggpmisc)

# # prok for breed ------------------------
# 
# source("/code/R_func/plot_PCoA.R")
# group <- read.delim("../profile/sample_group")
# profile <- readRDS("../profile/genomospecies.tpm.b50.rds")
# colors <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA")
# colors <- c(rep(colors, 3), "#E69F00")
# 
# group_x <- group %>%
#   filter(region == "Rectum" & group2 != "") %>%
#   select(sample, group = breed2, class) %>%
#   arrange(class, group) %>%
#   mutate(group = factor(group, unique(group)))
# 
# profile_x <- profile %>%
#   select(any_of(group_x$sample)) %>%
#   filter(rowSums(.) != 0)
# 
# # bray
# distance <- calcu_distance(profile_x, "bray")
# data <- calcu_PCoA(distance = distance, group = group_x, adonis2 = T)
# plot_data <- left_join(data$points, group_x, by = "sample")
# 
# ggscatterhist(plot_data, x = "PCoA1", y = "PCoA2", color = "group", shape = "class",
#               palette = colors, legend = "right", title = "Bray distance PCoA",
#               margin.plot = "boxplot", margin.plot.size = 3, main.plot.size = 1,
#               margin.params = list(outlier.shape = NA, fill = "group", width = 1),
#               xlab = data$eig_[1], ylab = data$eig_[2], subtitle = data$label)
# 
# ggsave("prok.beta_div.bray.for_breed.pdf", height = 7, width = 8)
# 
# plot_data <- reshape2::melt(as.matrix(distance)) %>%
#   filter(as.character(Var1) != as.character(Var2)) %>%
#   mutate_if(is.factor, as.character) %>%
#   left_join(group_x, by = c("Var1" = "sample"))
# 
# var <- plot_data %>%
#   filter(group == "JL-DLLW") %>%
#   pull(value) %>%
#   median()
# 
# ggerrorplot(plot_data, "group", "value", color = "class", xlab = "", x.text.angle = 90,
#             ylab = "Bray curtis-based dissimilarity", desc_stat = "median_q1q3",
#             legend = "none", palette = c("#56B4E9", "#E69F00"), ylim = c(0.5, NA)) +
#   stat_compare_means(ref.group = "JL-DLLW", label = "p.signif") +
#   geom_hline(yintercept = var, lty = "dashed", linewidth = .4)
# 
# ggsave("prok.beta_div.bray.for_breed.errorplot.pdf", height = 4, width = 6)
# 
# # jaccard
# distance <- calcu_distance(profile_x, "jaccard")
# data <- calcu_PCoA(distance = distance, group = group_x, adonis2 = T)
# plot_data <- left_join(data$points, group_x, by = "sample")
# 
# ggscatterhist(plot_data, x = "PCoA1", y = "PCoA2", color = "group", shape = "class",
#               palette = colors, legend = "right", title = "Jaccard distance PCoA",
#               margin.plot = "boxplot", margin.plot.size = 3, main.plot.size = 1,
#               margin.params = list(outlier.shape = NA, fill = "group", width = 1),
#               xlab = data$eig_[1], ylab = data$eig_[2], subtitle = data$label)
# 
# ggsave("prok.beta_div.jaccard.for_breed.pdf", height = 7, width = 8)
# 
# plot_data <- reshape2::melt(as.matrix(distance)) %>%
#   filter(as.character(Var1) != as.character(Var2)) %>%
#   mutate_if(is.factor, as.character) %>%
#   left_join(group_x, by = c("Var1" = "sample"))
# 
# var <- plot_data %>%
#   filter(group == "JL-DLLW") %>%
#   pull(value) %>%
#   median()
# 
# ggerrorplot(plot_data, "group", "value", color = "class", xlab = "", x.text.angle = 90,
#             ylab = "Jaccard-based dissimilarity", desc_stat = "median_q1q3",
#             legend = "none", palette = c("#56B4E9", "#E69F00"), ylim = c(0.7, NA)) +
#   stat_compare_means(ref.group = "JL-DLLW", label = "p.signif") +
#   geom_hline(yintercept = var, lty = "dashed", linewidth = .4)
# 
# ggsave("prok.beta_div.jaccard.for_breed.errorplot.pdf", height = 4, width = 6)
# 
# # weighted-unifrac
# tr <- ape::read.tree("../profile/genomospecies.tre")
# distance <- calcu_distance(profile_x, "unifrac", tr, weighted = T)
# data <- calcu_PCoA(distance = distance, group = group_x, adonis2 = T)
# plot_data <- left_join(data$points, group_x, by = "sample")
# 
# ggscatterhist(plot_data, x = "PCoA1", y = "PCoA2", color = "group", shape = "class",
#               palette = colors, legend = "right", title = "Weighted-unifrac distance PCoA",
#               margin.plot = "boxplot", margin.plot.size = 3, main.plot.size = 1,
#               margin.params = list(outlier.shape = NA, fill = "group", width = 1),
#               xlab = data$eig_[1], ylab = data$eig_[2], subtitle = data$label)
# 
# ggsave("prok.beta_div.weighted_unifrac.for_breed.pdf", height = 7, width = 8)
# 
# plot_data <- reshape2::melt(as.matrix(distance)) %>%
#   filter(as.character(Var1) != as.character(Var2)) %>%
#   mutate_if(is.factor, as.character) %>%
#   left_join(group_x, by = c("Var1" = "sample"))
# 
# var <- plot_data %>%
#   filter(group == "JL-DLLW") %>%
#   pull(value) %>%
#   median()
# 
# ggerrorplot(plot_data, "group", "value", color = "class", xlab = "", x.text.angle = 90,
#             ylab = "Weighted-unifrac-based dissimilarity", desc_stat = "median_q1q3",
#             legend = "none", palette = c("#56B4E9", "#E69F00"), ylim = c(NA, 0.43)) +
#   stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", label.y = 0.43) +
#   geom_hline(yintercept = var, lty = "dashed", linewidth = .4)
# 
# ggsave("prok.beta_div.weighted_unifrac.for_breed.errorplot.pdf", height = 4, width = 6)
# 
# # unweighted-unifrac
# distance <- calcu_distance(profile_x, "unifrac", tr, weighted = F)
# data <- calcu_PCoA(distance = distance, group = group_x, adonis2 = T)
# plot_data <- left_join(data$points, group_x, by = "sample")
# 
# ggscatterhist(plot_data, x = "PCoA1", y = "PCoA2", color = "group", shape = "class",
#               palette = colors, legend = "right", title = "Unweighted-unifrac distance PCoA",
#               margin.plot = "boxplot", margin.plot.size = 3, main.plot.size = 1,
#               margin.params = list(outlier.shape = NA, fill = "group", width = 1),
#               xlab = data$eig_[1], ylab = data$eig_[2], subtitle = data$label)
# 
# ggsave("prok.beta_div.unweighted_unifrac.for_breed.pdf", height = 7, width = 8)
# 
# plot_data <- reshape2::melt(as.matrix(distance)) %>%
#   filter(as.character(Var1) != as.character(Var2)) %>%
#   mutate_if(is.factor, as.character) %>%
#   left_join(group_x, by = c("Var1" = "sample"))
# 
# var <- plot_data %>%
#   filter(group == "JL-DLLW") %>%
#   pull(value) %>%
#   median()
# 
# ggerrorplot(plot_data, "group", "value", color = "class", xlab = "", x.text.angle = 90,
#             ylab = "Unweighted-unifrac-based dissimilarity", desc_stat = "median_q1q3",
#             legend = "none", palette = c("#56B4E9", "#E69F00"), ylim = c(0.5, NA)) +
#   stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", label.y = 0.8) +
#   geom_hline(yintercept = var, lty = "dashed", linewidth = .4)
# 
# ggsave("prok.beta_div.unweighted_unifrac.for_breed.errorplot.pdf", height = 4, width = 6)
# 
# # prok for region facet by region ------------------------------------------------------------------------------------
# 
# source("/code/R_func/plot_PCoA.R")
# colors_region <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#ffed6f","#7fc97f")
# colors_breed <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f")
# region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")
# 
# group_x <- group %>%
#   filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
#   select(sample, breed = breed2, region)
# 
# profile_x <- profile %>%
#   select(any_of(group_x$sample)) %>%
#   filter(rowSums(.) != 0)
# 
# # bray
# plot_list <- list()
# plot_list2 <- list()
# adonis2_data <- data.frame(name = character(), r2 = numeric(), pval = numeric())
# for (i in region_order) {
#   group_i <- filter(group_x, region == i) %>%
#     rename(group = breed)
#   profile_i <- select(profile, any_of(group_i$sample))
#   data <- calcu_PCoA(profile_i, group = group_i, adonis2 = T)
# 
#   plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line",
#                               title = paste0("Bray-distance PCoA | ", i),
#                               add_lab_to_plot = T, dis_method = "bray",
#                               group_color = colors_breed, show_legend = F)
# 
#   adonis2_data <- adonis2_data %>%
#     add_row(name = i, r2 = data$adonis2_r2, pval = data$adonis2_p)
# 
#   plot_data <- reshape2::melt(as.matrix(data$dist)) %>%
#     filter(as.character(Var1) != as.character(Var2)) %>%
#     mutate_if(is.factor, as.character) %>%
#     left_join(group_i, by = c("Var1" = "sample"))
# 
#   var <- plot_data %>%
#     filter(group == "JL-DLLW") %>%
#     pull(value) %>%
#     median()
# 
#   plot_list2[[i]] <-  ggerrorplot(plot_data, "group", "value", color = "group", xlab = "", x.text.angle = 90,
#               ylab = "Bray curtis-based dissimilarity", desc_stat = "median_q1q3",
#               legend = "none", palette = colors_region, title = paste0("Bray | ", i)) +
#     stat_compare_means(ref.group = "JL-DLLW", label = "p.signif") +
#     geom_hline(yintercept = var, lty = "dashed", linewidth = .4) +
#     theme(aspect.ratio = 1)
# }
# 
# cowplot::plot_grid(plotlist = plot_list, nrow = 1)
# ggsave("prok.beta_div.bray.for_region.facet_by_region.pdf", height = 3, width = 22)
# 
# cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
# ggsave("prok.beta_div.bray.for_region.facet_by_region.errorplot.pdf", height = 4, width = 20)
# 
# write.table(adonis2_data, "prok.beta_div.bray.for_region.facet_by_region.adonis.txt", sep = "\t", quote = F, row.names = F)
# 
# adonis2_data %>%
#   mutate(name = factor(name, region_order),
#          plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>%
#   ggbarplot(x = "name", y = "r2", label = .$plab, fill = "name", palette = colors_region,
#             legend = "none", color = NA, x.text.angle = 90)
# ggsave("prok.beta_div.bray.for_region.facet_by_region.adonis.barplot.pdf", height = 5, width = 5)
# 
# # jaccard
# plot_list <- list()
# plot_list2 <- list()
# adonis2_data <- data.frame(name = character(), r2 = numeric(), pval = numeric())
# for (i in region_order) {
#   group_i <- filter(group_x, region == i) %>%
#     rename(group = breed)
#   profile_i <- select(profile, any_of(group_i$sample))
#   data <- calcu_PCoA(profile_i, group = group_i, adonis2 = T, dis_method = "jaccard")
# 
#   plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line",
#                               title = paste0("Jaccard-distance PCoA | ", i),
#                               add_lab_to_plot = T, dis_method = "bray",
#                               group_color = colors_breed, show_legend = F)
# 
#   adonis2_data <- adonis2_data %>%
#     add_row(name = i, r2 = data$adonis2_r2, pval = data$adonis2_p)
# 
#   plot_data <- reshape2::melt(as.matrix(data$dist)) %>%
#     filter(as.character(Var1) != as.character(Var2)) %>%
#     mutate_if(is.factor, as.character) %>%
#     left_join(group_i, by = c("Var1" = "sample"))
# 
#   var <- plot_data %>%
#     filter(group == "JL-DLLW") %>%
#     pull(value) %>%
#     median()
# 
#   plot_list2[[i]] <-  ggerrorplot(plot_data, "group", "value", color = "group", xlab = "", x.text.angle = 90,
#                                   ylab = "Jaccard-based dissimilarity", desc_stat = "median_q1q3",
#                                   legend = "none", palette = colors_region, title = paste0("Jaccard | ", i)) +
#     stat_compare_means(ref.group = "JL-DLLW", label = "p.signif") +
#     geom_hline(yintercept = var, lty = "dashed", linewidth = .4) +
#     theme(aspect.ratio = 1)
# }
# 
# cowplot::plot_grid(plotlist = plot_list, nrow = 1)
# ggsave("prok.beta_div.jaccard.for_region.facet_by_region.pdf", height = 3, width = 22)
# 
# cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
# ggsave("prok.beta_div.jaccard.for_region.facet_by_region.errorplot.pdf", height = 4, width = 20)
# 
# write.table(adonis2_data, "prok.beta_div.jaccard.for_region.facet_by_region.adonis.txt", sep = "\t", quote = F, row.names = F)
# adonis2_data %>%
#   mutate(name = factor(name, region_order),
#          plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>%
#   ggbarplot(x = "name", y = "r2", label = .$plab, fill = "name", palette = colors_region,
#             legend = "none", color = NA, x.text.angle = 90)
# ggsave("prok.beta_div.jaccard.for_region.facet_by_region.adonis.barplot.pdf", height = 5, width = 5)
# 
# # weighted-unifrac
# tr <- ape::read.tree("../profile/genomospecies.tre")
# plot_list <- list()
# plot_list2 <- list()
# adonis2_data <- data.frame(name = character(), r2 = numeric(), pval = numeric())
# for (i in region_order) {
#   group_i <- filter(group_x, region == i) %>%
#     rename(group = breed)
#   profile_i <- select(profile, any_of(group_i$sample))
#   distance <- calcu_distance(profile_i, "unifrac", tr, weighted = T)
#   data <- calcu_PCoA(distance = distance, group = group_i, adonis2 = T)
# 
#   plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line",
#                               title = paste0("Weighted unifrac-distance PCoA | ", i),
#                               add_lab_to_plot = T, group_color = colors_breed, show_legend = F)
# 
#   adonis2_data <- adonis2_data %>%
#     add_row(name = i, r2 = data$adonis2_r2, pval = data$adonis2_p)
# 
#   plot_data <- reshape2::melt(as.matrix(data$dist)) %>%
#     filter(as.character(Var1) != as.character(Var2)) %>%
#     mutate_if(is.factor, as.character) %>%
#     left_join(group_i, by = c("Var1" = "sample"))
# 
#   var <- plot_data %>%
#     filter(group == "JL-DLLW") %>%
#     pull(value) %>%
#     median()
# 
#   plot_list2[[i]] <-  ggerrorplot(plot_data, "group", "value", color = "group", xlab = "", x.text.angle = 90,
#                                   ylab = "Weighted unifrac-based dissimilarity", desc_stat = "median_q1q3",
#                                   legend = "none", palette = colors_region, title = paste0("Weighted unifrac | ", i)) +
#     stat_compare_means(ref.group = "JL-DLLW", label = "p.signif") +
#     geom_hline(yintercept = var, lty = "dashed", linewidth = .4) +
#     theme(aspect.ratio = 1)
# }
# 
# cowplot::plot_grid(plotlist = plot_list, nrow = 1)
# ggsave("prok.beta_div.weighted_unifrac.for_region.facet_by_region.pdf", height = 3, width = 22)
# 
# cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
# ggsave("prok.beta_div.weighted_unifrac.for_region.facet_by_region.errorplot.pdf", height = 4, width = 20)
# 
# write.table(adonis2_data, "prok.beta_div.weighted_unifrac.for_region.facet_by_region.adonis.txt", sep = "\t", quote = F, row.names = F)
# adonis2_data %>%
#   mutate(name = factor(name, region_order),
#          plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>%
#   ggbarplot(x = "name", y = "r2", label = .$plab, fill = "name", palette = colors_region,
#             legend = "none", color = NA, x.text.angle = 90)
# ggsave("prok.beta_div.weighted_unifrac.for_region.facet_by_region.adonis.barplot.pdf", height = 5, width = 5)
# 
# # unweighted-unifrac
# plot_list <- list()
# plot_list2 <- list()
# adonis2_data <- data.frame(name = character(), r2 = numeric(), pval = numeric())
# for (i in region_order) {
#   group_i <- filter(group_x, region == i) %>%
#     rename(group = breed)
#   profile_i <- select(profile, any_of(group_i$sample))
#   distance <- calcu_distance(profile_i, "unifrac", tr, weighted = F)
#   data <- calcu_PCoA(distance = distance, group = group_i, adonis2 = T)
# 
#   plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line",
#                               title = paste0("Unweighted unifrac-distance PCoA | ", i),
#                               add_lab_to_plot = T, group_color = colors_breed, show_legend = F)
# 
#   adonis2_data <- adonis2_data %>%
#     add_row(name = i, r2 = data$adonis2_r2, pval = data$adonis2_p)
# 
#   plot_data <- reshape2::melt(as.matrix(data$dist)) %>%
#     filter(as.character(Var1) != as.character(Var2)) %>%
#     mutate_if(is.factor, as.character) %>%
#     left_join(group_i, by = c("Var1" = "sample"))
# 
#   var <- plot_data %>%
#     filter(group == "JL-DLLW") %>%
#     pull(value) %>%
#     median()
# 
#   plot_list2[[i]] <-  ggerrorplot(plot_data, "group", "value", color = "group", xlab = "", x.text.angle = 90,
#                                   ylab = "Unweighted unifrac-based dissimilarity", desc_stat = "median_q1q3",
#                                   legend = "none", palette = colors_region, title = paste0("Unweighted unifrac | ", i)) +
#     stat_compare_means(ref.group = "JL-DLLW", label = "p.signif") +
#     geom_hline(yintercept = var, lty = "dashed", linewidth = .4) +
#     theme(aspect.ratio = 1)
# }
# 
# cowplot::plot_grid(plotlist = plot_list, nrow = 1)
# ggsave("prok.beta_div.unweighted_unifrac.for_region.facet_by_region.pdf", height = 3, width = 22)
# 
# cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
# ggsave("prok.beta_div.unweighted_unifrac.for_region.facet_by_region.errorplot.pdf", height = 4, width = 20)
# 
# write.table(adonis2_data, "prok.beta_div.unweighted_unifrac.for_region.facet_by_region.adonis.txt", sep = "\t", quote = F, row.names = F)
# adonis2_data %>%
#   mutate(name = factor(name, region_order),
#          plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>%
#   ggbarplot(x = "name", y = "r2", label = .$plab, fill = "name", palette = colors_region,
#             legend = "none", color = NA, x.text.angle = 90)
# ggsave("prok.beta_div.unweighted_unifrac.for_region.facet_by_region.adonis.barplot.pdf", height = 5, width = 5)
# 
# # prok for region facet by breed ------------------------------------------------------------------------------------
# 
# source("/code/R_func/plot_PCoA.R")
# colors_region <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#ffed6f","#7fc97f")
# colors_breed <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f")
# region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")
# breed_order <- c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")
# 
# group_x <- group %>%
#   filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
#   select(sample, breed = breed2, region) %>%
#   mutate(region = factor(region, region_order))
# 
# profile_x <- profile %>%
#   select(any_of(group_x$sample)) %>%
#   filter(rowSums(.) != 0)
# 
# # bray
# plot_list <- list()
# plot_list2 <- list()
# adonis2_data <- data.frame(name = character(), r2 = numeric(), pval = numeric())
# for (i in breed_order) {
#   group_i <- filter(group_x, breed == i) %>%
#     rename(group = region)
#   profile_i <- select(profile, any_of(group_i$sample))
#   data <- calcu_PCoA(profile_i, group = group_i, adonis2 = T, dis_method = "bray")
# 
#   plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line",
#                               title = paste0("Bray-distance PCoA | ", i),
#                               add_lab_to_plot = T, dis_method = "bray",
#                               group_color = colors_region, show_legend = F)
# 
#   adonis2_data <- adonis2_data %>%
#     add_row(name = i, r2 = data$adonis2_r2, pval = data$adonis2_p)
# 
#   plot_data <- reshape2::melt(as.matrix(data$dist)) %>%
#     filter(as.character(Var1) != as.character(Var2)) %>%
#     mutate_if(is.factor, as.character) %>%
#     left_join(group_i, by = c("Var1" = "sample"))
# 
#   plot_list2[[i]] <- ggerrorplot(plot_data, "group", "value", color = "group", xlab = "", x.text.angle = 90,
#                                  ylab = "Bray curtis-based dissimilarity", desc_stat = "median_q1q3",
#                                  legend = "none", palette = colors_region, title = paste0("Bray | ", i)) +
#       geom_line(aes(group = 1, x = group, y = median),
#                 plot_data %>% group_by(group) %>% summarise(median = median(value)),
#                 lty = "dashed") +
#     theme(aspect.ratio = 2/3)
#   }
# 
# cowplot::plot_grid(plotlist = plot_list, nrow = 1)
# ggsave("prok.beta_div.bray.for_region.facet_by_breed.pdf", height = 3, width = 12)
# 
# cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
# ggsave("prok.beta_div.bray.for_region.facet_by_breed.errorplot.pdf", height = 4, width = 15)
# 
# write.table(adonis2_data, "prok.beta_div.bray.for_region.facet_by_breed.adonis.txt", sep = "\t", quote = F, row.names = F)
# adonis2_data %>%
#   mutate(name = factor(name, breed_order),
#          plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>%
#   ggbarplot(x = "name", y = "r2", label = .$plab, fill = "name", palette = colors_breed,
#             legend = "none", color = NA, x.text.angle = 90)
# ggsave("prok.beta_div.bray.for_region.facet_by_breed.adonis.barplot.pdf", height = 5, width = 5)
# 
# # jaccard
# plot_list <- list()
# plot_list2 <- list()
# adonis2_data <- data.frame(name = character(), r2 = numeric(), pval = numeric())
# for (i in breed_order) {
#   group_i <- filter(group_x, breed == i) %>%
#     rename(group = region)
#   profile_i <- select(profile, any_of(group_i$sample))
#   data <- calcu_PCoA(profile_i, group = group_i, adonis2 = T, dis_method = "jaccard")
# 
#   plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line",
#                               title = paste0("Jaccard-distance PCoA | ", i),
#                               add_lab_to_plot = T, dis_method = "jaccard",
#                               group_color = colors_region, show_legend = F)
# 
#   adonis2_data <- adonis2_data %>%
#     add_row(name = i, r2 = data$adonis2_r2, pval = data$adonis2_p)
# 
#   plot_data <- reshape2::melt(as.matrix(data$dist)) %>%
#     filter(as.character(Var1) != as.character(Var2)) %>%
#     mutate_if(is.factor, as.character) %>%
#     left_join(group_i, by = c("Var1" = "sample"))
# 
#   plot_list2[[i]] <- ggerrorplot(plot_data, "group", "value", color = "group", xlab = "", x.text.angle = 90,
#                                  ylab = "Jaccard-based dissimilarity", desc_stat = "median_q1q3",
#                                  legend = "none", palette = colors_region, title = paste0("Jaccard | ", i)) +
#     geom_line(aes(group = 1, x = group, y = median),
#               plot_data %>% group_by(group) %>% summarise(median = median(value)),
#               lty = "dashed") +
#     theme(aspect.ratio = 2/3)
# }
# 
# cowplot::plot_grid(plotlist = plot_list, nrow = 1)
# ggsave("prok.beta_div.jaccard.for_region.facet_by_breed.pdf", height = 3, width = 12)
# 
# cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
# ggsave("prok.beta_div.jaccard.for_region.facet_by_breed.errorplot.pdf", height = 4, width = 15)
# 
# write.table(adonis2_data, "prok.beta_div.jaccard.for_region.facet_by_breed.adonis.txt", sep = "\t", quote = F, row.names = F)
# adonis2_data %>%
#   mutate(name = factor(name, breed_order),
#          plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>%
#   ggbarplot(x = "name", y = "r2", label = .$plab, fill = "name", palette = colors_breed,
#             legend = "none", color = NA, x.text.angle = 90)
# ggsave("prok.beta_div.jaccard.for_region.facet_by_breed.adonis.barplot.pdf", height = 5, width = 5)
# 
# # weighted-unifrac
# tr <- ape::read.tree("../profile/genomospecies.tre")
# plot_list <- list()
# plot_list2 <- list()
# adonis2_data <- data.frame(name = character(), r2 = numeric(), pval = numeric())
# for (i in breed_order) {
#   group_i <- filter(group_x, breed == i) %>%
#     rename(group = region)
#   profile_i <- select(profile, any_of(group_i$sample))
#   distance <- calcu_distance(profile_i, "unifrac", tr, weighted = T)
#   data <- calcu_PCoA(distance = distance, group = group_i, adonis2 = T)
# 
#   plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line",
#                               title = paste0("Weighted-unifrac-distance PCoA | ", i),
#                               add_lab_to_plot = T, group_color = colors_region, show_legend = F)
# 
#   adonis2_data <- adonis2_data %>%
#     add_row(name = i, r2 = data$adonis2_r2, pval = data$adonis2_p)
# 
#   plot_data <- reshape2::melt(as.matrix(data$dist)) %>%
#     filter(as.character(Var1) != as.character(Var2)) %>%
#     mutate_if(is.factor, as.character) %>%
#     left_join(group_i, by = c("Var1" = "sample"))
# 
#   plot_list2[[i]] <- ggerrorplot(plot_data, "group", "value", color = "group", xlab = "", x.text.angle = 90,
#                                  ylab = "Weighted-unifrac-based dissimilarity", desc_stat = "median_q1q3",
#                                  legend = "none", palette = colors_region, title = paste0("Weighted-unifrac | ", i)) +
#     geom_line(aes(group = 1, x = group, y = median),
#               plot_data %>% group_by(group) %>% summarise(median = median(value)),
#               lty = "dashed") +
#     theme(aspect.ratio = 2/3)
# }
# 
# cowplot::plot_grid(plotlist = plot_list, nrow = 1)
# ggsave("prok.beta_div.weighted_unifrac.for_region.facet_by_breed.pdf", height = 3, width = 12)
# 
# cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
# ggsave("prok.beta_div.weighted_unifrac.for_region.facet_by_breed.errorplot.pdf", height = 4, width = 15)
# 
# write.table(adonis2_data, "prok.beta_div.weighted_unifrac.for_region.facet_by_breed.adonis.txt", sep = "\t", quote = F, row.names = F)
# adonis2_data %>%
#   mutate(name = factor(name, breed_order),
#          plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>%
#   ggbarplot(x = "name", y = "r2", label = .$plab, fill = "name", palette = colors_breed,
#             legend = "none", color = NA, x.text.angle = 90)
# ggsave("prok.beta_div.weighted_unifrac.for_region.facet_by_breed.adonis.barplot.pdf", height = 5, width = 5)
# 
# # unweighted-unifrac
# plot_list <- list()
# plot_list2 <- list()
# adonis2_data <- data.frame(name = character(), r2 = numeric(), pval = numeric())
# for (i in breed_order) {
#   group_i <- filter(group_x, breed == i) %>%
#     rename(group = region)
#   profile_i <- select(profile, any_of(group_i$sample))
#   distance <- calcu_distance(profile_i, "unifrac", tr, weighted = F)
#   data <- calcu_PCoA(distance = distance, group = group_i, adonis2 = T)
# 
#   plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line",
#                               title = paste0("Unweighted-unifrac-distance PCoA | ", i),
#                               add_lab_to_plot = T, group_color = colors_region, show_legend = F)
# 
#   adonis2_data <- adonis2_data %>%
#     add_row(name = i, r2 = data$adonis2_r2, pval = data$adonis2_p)
# 
#   plot_data <- reshape2::melt(as.matrix(data$dist)) %>%
#     filter(as.character(Var1) != as.character(Var2)) %>%
#     mutate_if(is.factor, as.character) %>%
#     left_join(group_i, by = c("Var1" = "sample"))
# 
#   plot_list2[[i]] <- ggerrorplot(plot_data, "group", "value", color = "group", xlab = "", x.text.angle = 90,
#                                  ylab = "Unweighted-unifrac-based dissimilarity", desc_stat = "median_q1q3",
#                                  legend = "none", palette = colors_region, title = paste0("Unweighted-unifrac | ", i)) +
#     geom_line(aes(group = 1, x = group, y = median),
#               plot_data %>% group_by(group) %>% summarise(median = median(value)),
#               lty = "dashed") +
#     theme(aspect.ratio = 2/3)
# }
# 
# cowplot::plot_grid(plotlist = plot_list, nrow = 1)
# ggsave("prok.beta_div.unweighted_unifrac.for_region.facet_by_breed.pdf", height = 3, width = 12)
# 
# cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
# ggsave("prok.beta_div.unweighted_unifrac.for_region.facet_by_breed.errorplot.pdf", height = 4, width = 15)
# 
# write.table(adonis2_data, "prok.beta_div.unweighted_unifrac.for_region.facet_by_breed.adonis.txt", sep = "\t", quote = F, row.names = F)
# adonis2_data %>%
#   mutate(name = factor(name, breed_order),
#          plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>%
#   ggbarplot(x = "name", y = "r2", label = .$plab, fill = "name", palette = colors_breed,
#             legend = "none", color = NA, x.text.angle = 90)
# ggsave("prok.beta_div.unweighted_unifrac.for_region.facet_by_breed.adonis.barplot.pdf", height = 5, width = 5)
# 
# # phage for breed ---------------------------------------------------------
# 
# source("/code/R_func/plot_PCoA.R")
# group <- read.delim("../profile/sample_group")
# profile <- readRDS("../profile/final_vOTUs.tpm.rds")
# colors <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA")
# colors <- c(rep(colors, 3), "#E69F00")
# 
# group_x <- group %>%
#   filter(region == "Rectum" & group2 != "") %>%
#   select(sample, group = breed2, class) %>%
#   arrange(class, group) %>%
#   mutate(group = factor(group, unique(group)))
# 
# profile_x <- profile %>%
#   select(any_of(group_x$sample)) %>%
#   filter(rowSums(.) != 0)
# 
# # bray
# distance <- calcu_distance(profile_x, "bray")
# data <- calcu_PCoA(distance = distance, group = group_x, adonis2 = T)
# plot_data <- left_join(data$points, group_x, by = "sample")
# 
# ggscatterhist(plot_data, x = "PCoA1", y = "PCoA2", color = "group", shape = "class",
#               palette = colors, legend = "right", title = "Bray distance PCoA",
#               margin.plot = "boxplot", margin.plot.size = 3, main.plot.size = 1,
#               margin.params = list(outlier.shape = NA, fill = "group", width = 1),
#               xlab = data$eig_[1], ylab = data$eig_[2], subtitle = data$label)
# 
# ggsave("phage.beta_div.bray.for_breed.pdf", height = 7, width = 8)
# 
# plot_data <- reshape2::melt(as.matrix(distance)) %>%
#   filter(as.character(Var1) != as.character(Var2)) %>%
#   mutate_if(is.factor, as.character) %>%
#   left_join(group_x, by = c("Var1" = "sample"))
# 
# var <- plot_data %>%
#   filter(group == "JL-DLLW") %>%
#   pull(value) %>%
#   median()
# 
# ggerrorplot(plot_data, "group", "value", color = "class", xlab = "", x.text.angle = 90,
#             ylab = "Bray curtis-based dissimilarity", desc_stat = "median_q1q3",
#             legend = "none", palette = c("#56B4E9", "#E69F00"), ylim = c(0.6, NA)) +
#   stat_compare_means(ref.group = "JL-DLLW", label = "p.signif") +
#   geom_hline(yintercept = var, lty = "dashed", linewidth = .4)
# 
# ggsave("phage.beta_div.bray.for_breed.errorplot.pdf", height = 4, width = 6)
# 
# # jaccard
# distance <- calcu_distance(profile_x, "jaccard")
# data <- calcu_PCoA(distance = distance, group = group_x, adonis2 = T)
# plot_data <- left_join(data$points, group_x, by = "sample")
# 
# ggscatterhist(plot_data, x = "PCoA1", y = "PCoA2", color = "group", shape = "class",
#               palette = colors, legend = "right", title = "Jaccard distance PCoA",
#               margin.plot = "boxplot", margin.plot.size = 3, main.plot.size = 1,
#               margin.params = list(outlier.shape = NA, fill = "group", width = 1),
#               xlab = data$eig_[1], ylab = data$eig_[2], subtitle = data$label)
# 
# ggsave("phage.beta_div.jaccard.for_breed.pdf", height = 7, width = 8)
# 
# plot_data <- reshape2::melt(as.matrix(distance)) %>%
#   filter(as.character(Var1) != as.character(Var2)) %>%
#   mutate_if(is.factor, as.character) %>%
#   left_join(group_x, by = c("Var1" = "sample"))
# 
# var <- plot_data %>%
#   filter(group == "JL-DLLW") %>%
#   pull(value) %>%
#   median()
# 
# ggerrorplot(plot_data, "group", "value", color = "class", xlab = "", x.text.angle = 90,
#             ylab = "Jaccard-based dissimilarity", desc_stat = "median_q1q3",
#             legend = "none", palette = c("#56B4E9", "#E69F00"), ylim = c(0.79, NA)) +
#   stat_compare_means(ref.group = "JL-DLLW", label = "p.signif") +
#   geom_hline(yintercept = var, lty = "dashed", linewidth = .4)
# 
# ggsave("phage.beta_div.jaccard.for_breed.errorplot.pdf", height = 4, width = 6)
# 
# # weighted-unifrac
# tr <- ape::read.tree("../profile/genomospecies.tre")
# distance <- calcu_distance(profile_x, "unifrac", tr, weighted = T)
# data <- calcu_PCoA(distance = distance, group = group_x, adonis2 = T)
# plot_data <- left_join(data$points, group_x, by = "sample")
# 
# ggscatterhist(plot_data, x = "PCoA1", y = "PCoA2", color = "group", shape = "class",
#               palette = colors, legend = "right", title = "Weighted-unifrac distance PCoA",
#               margin.plot = "boxplot", margin.plot.size = 3, main.plot.size = 1,
#               margin.params = list(outlier.shape = NA, fill = "group", width = 1),
#               xlab = data$eig_[1], ylab = data$eig_[2], subtitle = data$label)
# 
# ggsave("phage.beta_div.weighted_unifrac.for_breed.pdf", height = 7, width = 8)
# 
# plot_data <- reshape2::melt(as.matrix(distance)) %>%
#   filter(as.character(Var1) != as.character(Var2)) %>%
#   mutate_if(is.factor, as.character) %>%
#   left_join(group_x, by = c("Var1" = "sample"))
# 
# var <- plot_data %>%
#   filter(group == "JL-DLLW") %>%
#   pull(value) %>%
#   median()
# 
# ggerrorplot(plot_data, "group", "value", color = "class", xlab = "", x.text.angle = 90,
#             ylab = "Weighted-unifrac-based dissimilarity", desc_stat = "median_q1q3",
#             legend = "none", palette = c("#56B4E9", "#E69F00"), ylim = c(NA, 0.43)) +
#   stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", label.y = 0.43) +
#   geom_hline(yintercept = var, lty = "dashed", linewidth = .4)
# 
# ggsave("phage.beta_div.weighted_unifrac.for_breed.errorplot.pdf", height = 4, width = 6)
# 
# # unweighted-unifrac
# distance <- calcu_distance(profile_x, "unifrac", tr, weighted = F)
# data <- calcu_PCoA(distance = distance, group = group_x, adonis2 = T)
# plot_data <- left_join(data$points, group_x, by = "sample")
# 
# ggscatterhist(plot_data, x = "PCoA1", y = "PCoA2", color = "group", shape = "class",
#               palette = colors, legend = "right", title = "Unweighted-unifrac distance PCoA",
#               margin.plot = "boxplot", margin.plot.size = 3, main.plot.size = 1,
#               margin.params = list(outlier.shape = NA, fill = "group", width = 1),
#               xlab = data$eig_[1], ylab = data$eig_[2], subtitle = data$label)
# 
# ggsave("phage.beta_div.unweighted_unifrac.for_breed.pdf", height = 7, width = 8)
# 
# plot_data <- reshape2::melt(as.matrix(distance)) %>%
#   filter(as.character(Var1) != as.character(Var2)) %>%
#   mutate_if(is.factor, as.character) %>%
#   left_join(group_x, by = c("Var1" = "sample"))
# 
# var <- plot_data %>%
#   filter(group == "JL-DLLW") %>%
#   pull(value) %>%
#   median()
# 
# ggerrorplot(plot_data, "group", "value", color = "class", xlab = "", x.text.angle = 90,
#             ylab = "Unweighted-unifrac-based dissimilarity", desc_stat = "median_q1q3",
#             legend = "none", palette = c("#56B4E9", "#E69F00"), ylim = c(0.5, NA)) +
#   stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", label.y = 0.8) +
#   geom_hline(yintercept = var, lty = "dashed", linewidth = .4)
# 
# ggsave("phage.beta_div.unweighted_unifrac.for_breed.errorplot.pdf", height = 4, width = 6)
# 
# # phage for region facet by region ------------------------------------------------------------------------------------
# 
# source("/code/R_func/plot_PCoA.R")
# colors_region <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#ffed6f","#7fc97f")
# colors_breed <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f")
# region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")
# 
# group_x <- group %>%
#   filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
#   select(sample, breed = breed2, region)
# 
# profile_x <- profile %>%
#   select(any_of(group_x$sample)) %>%
#   filter(rowSums(.) != 0)
# 
# # bray
# plot_list <- list()
# plot_list2 <- list()
# adonis2_data <- data.frame(name = character(), r2 = numeric(), pval = numeric())
# for (i in region_order) {
#   group_i <- filter(group_x, region == i) %>%
#     rename(group = breed)
#   profile_i <- select(profile, any_of(group_i$sample))
#   data <- calcu_PCoA(profile_i, group = group_i, adonis2 = T)
# 
#   plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line",
#                               title = paste0("Bray-distance PCoA | ", i),
#                               add_lab_to_plot = T, dis_method = "bray",
#                               group_color = colors_breed, show_legend = F)
# 
#   adonis2_data <- adonis2_data %>%
#     add_row(name = i, r2 = data$adonis2_r2, pval = data$adonis2_p)
# 
#   plot_data <- reshape2::melt(as.matrix(data$dist)) %>%
#     filter(as.character(Var1) != as.character(Var2)) %>%
#     mutate_if(is.factor, as.character) %>%
#     left_join(group_i, by = c("Var1" = "sample"))
# 
#   var <- plot_data %>%
#     filter(group == "JL-DLLW") %>%
#     pull(value) %>%
#     median()
# 
#   plot_list2[[i]] <-  ggerrorplot(plot_data, "group", "value", color = "group", xlab = "", x.text.angle = 90,
#                                   ylab = "Bray curtis-based dissimilarity", desc_stat = "median_q1q3",
#                                   legend = "none", palette = colors_region, title = paste0("Bray | ", i)) +
#     stat_compare_means(ref.group = "JL-DLLW", label = "p.signif") +
#     geom_hline(yintercept = var, lty = "dashed", linewidth = .4) +
#     theme(aspect.ratio = 1)
# }
# 
# cowplot::plot_grid(plotlist = plot_list, nrow = 1)
# ggsave("phage.beta_div.bray.for_region.facet_by_region.pdf", height = 3, width = 22)
# 
# cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
# ggsave("phage.beta_div.bray.for_region.facet_by_region.errorplot.pdf", height = 4, width = 20)
# 
# write.table(adonis2_data, "phage.beta_div.bray.for_region.facet_by_region.adonis.txt", sep = "\t", quote = F, row.names = F)
# adonis2_data %>%
#   mutate(name = factor(name, region_order),
#          plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>%
#   ggbarplot(x = "name", y = "r2", label = .$plab, fill = "name", palette = colors_region,
#             legend = "none", color = NA, x.text.angle = 90)
# ggsave("phage.beta_div.bray.for_region.facet_by_region.adonis.barplot.pdf", height = 5, width = 5)
# 
# # jaccard
# plot_list <- list()
# plot_list2 <- list()
# adonis2_data <- data.frame(name = character(), r2 = numeric(), pval = numeric())
# for (i in region_order) {
#   group_i <- filter(group_x, region == i) %>%
#     rename(group = breed)
#   profile_i <- select(profile, any_of(group_i$sample))
#   data <- calcu_PCoA(profile_i, group = group_i, adonis2 = T, dis_method = "jaccard")
# 
#   plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line",
#                               title = paste0("Jaccard-distance PCoA | ", i),
#                               add_lab_to_plot = T, dis_method = "bray",
#                               group_color = colors_breed, show_legend = F)
# 
#   adonis2_data <- adonis2_data %>%
#     add_row(name = i, r2 = data$adonis2_r2, pval = data$adonis2_p)
# 
#   plot_data <- reshape2::melt(as.matrix(data$dist)) %>%
#     filter(as.character(Var1) != as.character(Var2)) %>%
#     mutate_if(is.factor, as.character) %>%
#     left_join(group_i, by = c("Var1" = "sample"))
# 
#   var <- plot_data %>%
#     filter(group == "JL-DLLW") %>%
#     pull(value) %>%
#     median()
# 
#   plot_list2[[i]] <-  ggerrorplot(plot_data, "group", "value", color = "group", xlab = "", x.text.angle = 90,
#                                   ylab = "Jaccard-based dissimilarity", desc_stat = "median_q1q3",
#                                   legend = "none", palette = colors_region, title = paste0("Jaccard | ", i)) +
#     stat_compare_means(ref.group = "JL-DLLW", label = "p.signif") +
#     geom_hline(yintercept = var, lty = "dashed", linewidth = .4) +
#     theme(aspect.ratio = 1)
# }
# 
# cowplot::plot_grid(plotlist = plot_list, nrow = 1)
# ggsave("phage.beta_div.jaccard.for_region.facet_by_region.pdf", height = 3, width = 22)
# 
# cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
# ggsave("phage.beta_div.jaccard.for_region.facet_by_region.errorplot.pdf", height = 4, width = 20)
# 
# write.table(adonis2_data, "phage.beta_div.jaccard.for_region.facet_by_region.adonis.txt", sep = "\t", quote = F, row.names = F)
# adonis2_data %>%
#   mutate(name = factor(name, region_order),
#          plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>%
#   ggbarplot(x = "name", y = "r2", label = .$plab, fill = "name", palette = colors_region,
#             legend = "none", color = NA, x.text.angle = 90)
# ggsave("phage.beta_div.jaccard.for_region.facet_by_region.adonis.barplot.pdf", height = 5, width = 5)
# 
# # weighted-unifrac
# tr <- ape::read.tree("../profile/genomospecies.tre")
# plot_list <- list()
# plot_list2 <- list()
# adonis2_data <- data.frame(name = character(), r2 = numeric(), pval = numeric())
# for (i in region_order) {
#   group_i <- filter(group_x, region == i) %>%
#     rename(group = breed)
#   profile_i <- select(profile, any_of(group_i$sample))
#   distance <- calcu_distance(profile_i, "unifrac", tr, weighted = T)
#   data <- calcu_PCoA(distance = distance, group = group_i, adonis2 = T)
# 
#   plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line",
#                               title = paste0("Weighted unifrac-distance PCoA | ", i),
#                               add_lab_to_plot = T, group_color = colors_breed, show_legend = F)
# 
#   adonis2_data <- adonis2_data %>%
#     add_row(name = i, r2 = data$adonis2_r2, pval = data$adonis2_p)
# 
#   plot_data <- reshape2::melt(as.matrix(data$dist)) %>%
#     filter(as.character(Var1) != as.character(Var2)) %>%
#     mutate_if(is.factor, as.character) %>%
#     left_join(group_i, by = c("Var1" = "sample"))
# 
#   var <- plot_data %>%
#     filter(group == "JL-DLLW") %>%
#     pull(value) %>%
#     median()
# 
#   plot_list2[[i]] <-  ggerrorplot(plot_data, "group", "value", color = "group", xlab = "", x.text.angle = 90,
#                                   ylab = "Weighted unifrac-based dissimilarity", desc_stat = "median_q1q3",
#                                   legend = "none", palette = colors_region, title = paste0("Weighted unifrac | ", i)) +
#     stat_compare_means(ref.group = "JL-DLLW", label = "p.signif") +
#     geom_hline(yintercept = var, lty = "dashed", linewidth = .4) +
#     theme(aspect.ratio = 1)
# }
# 
# cowplot::plot_grid(plotlist = plot_list, nrow = 1)
# ggsave("phage.beta_div.weighted_unifrac.for_region.facet_by_region.pdf", height = 3, width = 22)
# 
# cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
# ggsave("phage.beta_div.weighted_unifrac.for_region.facet_by_region.errorplot.pdf", height = 4, width = 20)
# 
# write.table(adonis2_data, "phage.beta_div.weighted_unifrac.for_region.facet_by_region.adonis.txt", sep = "\t", quote = F, row.names = F)
# adonis2_data %>%
#   mutate(name = factor(name, region_order),
#          plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>%
#   ggbarplot(x = "name", y = "r2", label = .$plab, fill = "name", palette = colors_region,
#             legend = "none", color = NA, x.text.angle = 90)
# ggsave("phage.beta_div.weighted_unifrac.for_region.facet_by_region.adonis.barplot.pdf", height = 5, width = 5)
# 
# # unweighted-unifrac
# plot_list <- list()
# plot_list2 <- list()
# adonis2_data <- data.frame(name = character(), r2 = numeric(), pval = numeric())
# for (i in region_order) {
#   group_i <- filter(group_x, region == i) %>%
#     rename(group = breed)
#   profile_i <- select(profile, any_of(group_i$sample))
#   distance <- calcu_distance(profile_i, "unifrac", tr, weighted = F)
#   data <- calcu_PCoA(distance = distance, group = group_i, adonis2 = T)
# 
#   plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line",
#                               title = paste0("Unweighted unifrac-distance PCoA | ", i),
#                               add_lab_to_plot = T, group_color = colors_breed, show_legend = F)
# 
#   adonis2_data <- adonis2_data %>%
#     add_row(name = i, r2 = data$adonis2_r2, pval = data$adonis2_p)
# 
#   plot_data <- reshape2::melt(as.matrix(data$dist)) %>%
#     filter(as.character(Var1) != as.character(Var2)) %>%
#     mutate_if(is.factor, as.character) %>%
#     left_join(group_i, by = c("Var1" = "sample"))
# 
#   var <- plot_data %>%
#     filter(group == "JL-DLLW") %>%
#     pull(value) %>%
#     median()
# 
#   plot_list2[[i]] <-  ggerrorplot(plot_data, "group", "value", color = "group", xlab = "", x.text.angle = 90,
#                                   ylab = "Unweighted unifrac-based dissimilarity", desc_stat = "median_q1q3",
#                                   legend = "none", palette = colors_region, title = paste0("Unweighted unifrac | ", i)) +
#     stat_compare_means(ref.group = "JL-DLLW", label = "p.signif") +
#     geom_hline(yintercept = var, lty = "dashed", linewidth = .4) +
#     theme(aspect.ratio = 1)
# }
# 
# cowplot::plot_grid(plotlist = plot_list, nrow = 1)
# ggsave("phage.beta_div.unweighted_unifrac.for_region.facet_by_region.pdf", height = 3, width = 22)
# 
# cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
# ggsave("phage.beta_div.unweighted_unifrac.for_region.facet_by_region.errorplot.pdf", height = 4, width = 20)
# 
# write.table(adonis2_data, "phage.beta_div.unweighted_unifrac.for_region.facet_by_region.adonis.txt", sep = "\t", quote = F, row.names = F)
# adonis2_data %>%
#   mutate(name = factor(name, region_order),
#          plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>%
#   ggbarplot(x = "name", y = "r2", label = .$plab, fill = "name", palette = colors_region,
#             legend = "none", color = NA, x.text.angle = 90)
# ggsave("phage.beta_div.unweighted_unifrac.for_region.facet_by_region.adonis.barplot.pdf", height = 5, width = 5)
# 
# # phage for region facet by breed ------------------------------------------------------------------------------------
# 
# source("/code/R_func/plot_PCoA.R")
# colors_region <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#ffed6f","#7fc97f")
# colors_breed <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f")
# region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")
# breed_order <- c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")
# 
# group_x <- group %>%
#   filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
#   select(sample, breed = breed2, region) %>%
#   mutate(region = factor(region, region_order))
# 
# profile_x <- profile %>%
#   select(any_of(group_x$sample)) %>%
#   filter(rowSums(.) != 0)
# 
# # bray
# plot_list <- list()
# plot_list2 <- list()
# adonis2_data <- data.frame(name = character(), r2 = numeric(), pval = numeric())
# for (i in breed_order) {
#   group_i <- filter(group_x, breed == i) %>%
#     rename(group = region)
#   profile_i <- select(profile, any_of(group_i$sample))
#   data <- calcu_PCoA(profile_i, group = group_i, adonis2 = T, dis_method = "bray")
# 
#   plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line",
#                               title = paste0("Bray-distance PCoA | ", i),
#                               add_lab_to_plot = T, dis_method = "bray",
#                               group_color = colors_region, show_legend = F)
# 
#   adonis2_data <- adonis2_data %>%
#     add_row(name = i, r2 = data$adonis2_r2, pval = data$adonis2_p)
# 
#   plot_data <- reshape2::melt(as.matrix(data$dist)) %>%
#     filter(as.character(Var1) != as.character(Var2)) %>%
#     mutate_if(is.factor, as.character) %>%
#     left_join(group_i, by = c("Var1" = "sample"))
# 
#   plot_list2[[i]] <- ggerrorplot(plot_data, "group", "value", color = "group", xlab = "", x.text.angle = 90,
#                                  ylab = "Bray curtis-based dissimilarity", desc_stat = "median_q1q3",
#                                  legend = "none", palette = colors_region, title = paste0("Bray | ", i)) +
#     geom_line(aes(group = 1, x = group, y = median),
#               plot_data %>% group_by(group) %>% summarise(median = median(value)),
#               lty = "dashed") +
#     theme(aspect.ratio = 2/3)
# }
# 
# cowplot::plot_grid(plotlist = plot_list, nrow = 1)
# ggsave("phage.beta_div.bray.for_region.facet_by_breed.pdf", height = 3, width = 12)
# 
# cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
# ggsave("phage.beta_div.bray.for_region.facet_by_breed.errorplot.pdf", height = 4, width = 15)
# 
# write.table(adonis2_data, "phage.beta_div.bray.for_region.facet_by_breed.adonis.txt", sep = "\t", quote = F, row.names = F)
# adonis2_data %>%
#   mutate(name = factor(name, breed_order),
#          plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>%
#   ggbarplot(x = "name", y = "r2", label = .$plab, fill = "name", palette = colors_breed,
#             legend = "none", color = NA, x.text.angle = 90)
# ggsave("phage.beta_div.bray.for_region.facet_by_breed.adonis.barplot.pdf", height = 5, width = 5)
# 
# # jaccard
# plot_list <- list()
# plot_list2 <- list()
# adonis2_data <- data.frame(name = character(), r2 = numeric(), pval = numeric())
# for (i in breed_order) {
#   group_i <- filter(group_x, breed == i) %>%
#     rename(group = region)
#   profile_i <- select(profile, any_of(group_i$sample))
#   data <- calcu_PCoA(profile_i, group = group_i, adonis2 = T, dis_method = "jaccard")
# 
#   plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line",
#                               title = paste0("Jaccard-distance PCoA | ", i),
#                               add_lab_to_plot = T, dis_method = "jaccard",
#                               group_color = colors_region, show_legend = F)
# 
#   adonis2_data <- adonis2_data %>%
#     add_row(name = i, r2 = data$adonis2_r2, pval = data$adonis2_p)
# 
#   plot_data <- reshape2::melt(as.matrix(data$dist)) %>%
#     filter(as.character(Var1) != as.character(Var2)) %>%
#     mutate_if(is.factor, as.character) %>%
#     left_join(group_i, by = c("Var1" = "sample"))
# 
#   plot_list2[[i]] <- ggerrorplot(plot_data, "group", "value", color = "group", xlab = "", x.text.angle = 90,
#                                  ylab = "Jaccard-based dissimilarity", desc_stat = "median_q1q3",
#                                  legend = "none", palette = colors_region, title = paste0("Jaccard | ", i)) +
#     geom_line(aes(group = 1, x = group, y = median),
#               plot_data %>% group_by(group) %>% summarise(median = median(value)),
#               lty = "dashed") +
#     theme(aspect.ratio = 2/3)
# }
# 
# cowplot::plot_grid(plotlist = plot_list, nrow = 1)
# ggsave("phage.beta_div.jaccard.for_region.facet_by_breed.pdf", height = 3, width = 12)
# 
# cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
# ggsave("phage.beta_div.jaccard.for_region.facet_by_breed.errorplot.pdf", height = 4, width = 15)
# 
# write.table(adonis2_data, "phage.beta_div.jaccard.for_region.facet_by_breed.adonis.txt", sep = "\t", quote = F, row.names = F)
# adonis2_data %>%
#   mutate(name = factor(name, breed_order),
#          plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>%
#   ggbarplot(x = "name", y = "r2", label = .$plab, fill = "name", palette = colors_breed,
#             legend = "none", color = NA, x.text.angle = 90)
# ggsave("phage.beta_div.jaccard.for_region.facet_by_breed.adonis.barplot.pdf", height = 5, width = 5)
# 
# # weighted-unifrac
# tr <- ape::read.tree("../profile/genomospecies.tre")
# plot_list <- list()
# plot_list2 <- list()
# adonis2_data <- data.frame(name = character(), r2 = numeric(), pval = numeric())
# for (i in breed_order) {
#   group_i <- filter(group_x, breed == i) %>%
#     rename(group = region)
#   profile_i <- select(profile, any_of(group_i$sample))
#   distance <- calcu_distance(profile_i, "unifrac", tr, weighted = T)
#   data <- calcu_PCoA(distance = distance, group = group_i, adonis2 = T)
# 
#   plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line",
#                               title = paste0("Weighted-unifrac-distance PCoA | ", i),
#                               add_lab_to_plot = T, group_color = colors_region, show_legend = F)
# 
#   adonis2_data <- adonis2_data %>%
#     add_row(name = i, r2 = data$adonis2_r2, pval = data$adonis2_p)
# 
#   plot_data <- reshape2::melt(as.matrix(data$dist)) %>%
#     filter(as.character(Var1) != as.character(Var2)) %>%
#     mutate_if(is.factor, as.character) %>%
#     left_join(group_i, by = c("Var1" = "sample"))
# 
#   plot_list2[[i]] <- ggerrorplot(plot_data, "group", "value", color = "group", xlab = "", x.text.angle = 90,
#                                  ylab = "Weighted-unifrac-based dissimilarity", desc_stat = "median_q1q3",
#                                  legend = "none", palette = colors_region, title = paste0("Weighted-unifrac | ", i)) +
#     geom_line(aes(group = 1, x = group, y = median),
#               plot_data %>% group_by(group) %>% summarise(median = median(value)),
#               lty = "dashed") +
#     theme(aspect.ratio = 2/3)
# }
# 
# cowplot::plot_grid(plotlist = plot_list, nrow = 1)
# ggsave("phage.beta_div.weighted_unifrac.for_region.facet_by_breed.pdf", height = 3, width = 12)
# 
# cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
# ggsave("phage.beta_div.weighted_unifrac.for_region.facet_by_breed.errorplot.pdf", height = 4, width = 15)
# 
# write.table(adonis2_data, "phage.beta_div.weighted_unifrac.for_region.facet_by_breed.adonis.txt", sep = "\t", quote = F, row.names = F)
# adonis2_data %>%
#   mutate(name = factor(name, breed_order),
#          plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>%
#   ggbarplot(x = "name", y = "r2", label = .$plab, fill = "name", palette = colors_breed,
#             legend = "none", color = NA, x.text.angle = 90)
# ggsave("phage.beta_div.weighted_unifrac.for_region.facet_by_breed.adonis.barplot.pdf", height = 5, width = 5)
# 
# # unweighted-unifrac
# plot_list <- list()
# plot_list2 <- list()
# adonis2_data <- data.frame(name = character(), r2 = numeric(), pval = numeric())
# for (i in breed_order) {
#   group_i <- filter(group_x, breed == i) %>%
#     rename(group = region)
#   profile_i <- select(profile, any_of(group_i$sample))
#   distance <- calcu_distance(profile_i, "unifrac", tr, weighted = F)
#   data <- calcu_PCoA(distance = distance, group = group_i, adonis2 = T)
# 
#   plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line",
#                               title = paste0("Unweighted-unifrac-distance PCoA | ", i),
#                               add_lab_to_plot = T, group_color = colors_region, show_legend = F)
# 
#   adonis2_data <- adonis2_data %>%
#     add_row(name = i, r2 = data$adonis2_r2, pval = data$adonis2_p)
# 
#   plot_data <- reshape2::melt(as.matrix(data$dist)) %>%
#     filter(as.character(Var1) != as.character(Var2)) %>%
#     mutate_if(is.factor, as.character) %>%
#     left_join(group_i, by = c("Var1" = "sample"))
# 
#   plot_list2[[i]] <- ggerrorplot(plot_data, "group", "value", color = "group", xlab = "", x.text.angle = 90,
#                                  ylab = "Unweighted-unifrac-based dissimilarity", desc_stat = "median_q1q3",
#                                  legend = "none", palette = colors_region, title = paste0("Unweighted-unifrac | ", i)) +
#     geom_line(aes(group = 1, x = group, y = median),
#               plot_data %>% group_by(group) %>% summarise(median = median(value)),
#               lty = "dashed") +
#     theme(aspect.ratio = 2/3)
# }
# 
# cowplot::plot_grid(plotlist = plot_list, nrow = 1)
# ggsave("phage.beta_div.unweighted_unifrac.for_region.facet_by_breed.pdf", height = 3, width = 12)
# 
# cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
# ggsave("phage.beta_div.unweighted_unifrac.for_region.facet_by_breed.errorplot.pdf", height = 4, width = 15)
# 
# write.table(adonis2_data, "phage.beta_div.unweighted_unifrac.for_region.facet_by_breed.adonis.txt", sep = "\t", quote = F, row.names = F)
# adonis2_data %>%
#   mutate(name = factor(name, breed_order),
#          plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>%
#   ggbarplot(x = "name", y = "r2", label = .$plab, fill = "name", palette = colors_breed,
#             legend = "none", color = NA, x.text.angle = 90)
# ggsave("phage.beta_div.unweighted_unifrac.for_region.facet_by_breed.adonis.barplot.pdf", height = 5, width = 5)
# 

# bray prok for breed ---------------------------------

library(patchwork)
source("/code/R_func/plot_PCoA.R")
source("/code/R_func/difference_analysis.R")
group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/genomospecies.tpm.b50.rds")
colors <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#F9C851","#ACDEF3","#f881bf",
            "#94d294","#d49264","#95beda","#eb9ecd","#F9C6B3","#fff087","#c3e6dd","#b3b3b3")

group_x <- group %>%
  filter(region == "Rectum" & group2 != "") %>% 
  select(sample, group = breed2, class) %>% 
  arrange(class, group) %>% 
  mutate(group = factor(group, unique(group)))

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

data <- calcu_PCoA(profile = profile_x, group = group_x, dis_method = "bray", adonis2 = T)

disper <- betadisper(data$dist, data$group$group) %>% 
  TukeyHSD() 
disper$group %>% 
  data.frame() %>% 
  rownames_to_column("group_pair") %>% 
  write.table("beta/prok.beta_div.bray.for_breed.disper.txt", sep = "\t", quote = F, row.names = F)

plot_data <- left_join(data$points, group_x, by = "sample")

p1 <- ggscatter(plot_data, x = "PCoA1", y = "PCoA2", color = "group", shape = "class",
          palette = colors, legend = "right", title = "Bray-Curtis distance-based PCoA",
          xlab = data$eig_[1], ylab = data$eig_[2], subtitle = data$label)

p2 <- ggboxplot(plot_data, x = "group", y = "PCoA1", fill = "group",
          palette = colors, rotate = T, outlier.shape = NA, color = "group",
          legend = "none") +
  stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", color = "red")

p3 <- ggboxplot(plot_data, x = "group", y = "PCoA2", fill = "group",
          palette = colors, x.text.angle = 90, outlier.shape = NA, color = "group",
          legend = "none") +
  stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", color = "red")

(p2 / p1) | (p2 / p3)

ggsave("beta/prok.beta_div.bray.for_breed.manual.pdf", height = 12, width = 14)

# adonis 分析
data <- calcu_pairwise_adonis(profile_x, group_x, dis_method = "bray")

write.table(data, "beta/prok.beta_div.bray.for_breed.adonis.tsv", sep = "\t", quote = F, row.names = F)

plot_data <- data %>% 
  filter(grepl("JL-DLLW", group_pair)) %>% 
  mutate(plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", ""))),
         name = gsub("_vs_JL-DLLW", "", group_pair),
         name = factor(name, name))

ggbarplot(plot_data, x = "name", y = "r2adj", fill = "name", palette = colors, legend = "none",
          x.text.angle = 90, xlab = "", ylab = "Adjusted R2", label = plot_data$plab, 
          lab.col = "red", width = .6)

ggsave("beta/prok.beta_div.bray.for_breed.adonis.pdf", height = 6, width = 6)

levels <- levels(group_x$group)

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

ggsave("beta/prok.beta_div.bray.for_breed.adonis2.pdf", width = 6, height = 6)

# bray prok for region facet by region -----------------

source("/code/R_func/plot_PCoA.R")
colors_region <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#F9C851","#7fc97f")
colors_breed <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f")
region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")

group_x <- group %>% 
  filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
  select(sample, breed = breed2, region)

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

# bray
plot_list <- list()
plot_list2 <- list()
adonis2_data <- data.frame(name = character(), r2 = numeric(), r2adj = numeric(), pval = numeric())
disper_data <- data.frame(name = character(), pval = numeric(), method = character())
disper_data2 <- rbind()

for (i in region_order) {
  group_i <- filter(group_x, region == i) %>% 
    rename(group = breed)
  profile_i <- select(profile, any_of(group_i$sample))
  data <- calcu_PCoA(profile_i, group = group_i, adonis2 = T)
  
  plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line", 
                              title = paste0("Bray-distance PCoA | ", i),
                              add_lab_to_plot = T, dis_method = "bray", 
                              group_color = colors_breed, show_legend = F)
  
  adonis2_data <- adonis2_data %>% 
    add_row(name = i, r2 = data$adonis2_r2, r2adj = data$adonis2_r2adj, pval = data$adonis2_p)
  
  disper <- betadisper(data$dist, data$group$group)
  
  disper_data <- disper_data %>% 
    add_row(name = i, pval = anova(disper)[1, 5], method = "anova")
  
  disper_data2 <- rbind(disper_data2, 
                        TukeyHSD(disper)[["group"]] %>% 
                          data.frame %>% 
                          rownames_to_column("group_pair") %>% 
                          add_column(group = i, method = "TukeyHSD"))

  plot_data <- reshape2::melt(as.matrix(data$dist)) %>% 
    filter(as.character(Var1) != as.character(Var2)) %>%
    mutate_if(is.factor, as.character) %>% 
    left_join(group_i, by = c("Var1" = "sample"))
  
  var <- plot_data %>% 
    filter(group == "JL-DLLW") %>% 
    pull(value) %>% 
    median()
  
  plot_list2[[i]] <- ggboxplot(plot_data, "group", "value", fill = "group", xlab = "", x.text.angle = 90,
                               ylab = "Bray curtis-based dissimilarity", 
                               legend = "none", palette = colors_region, title = paste0("Bray | ", i)) + 
    stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", color = "red") +
    geom_hline(yintercept = var, lty = "dashed", linewidth = .4) +
    theme(aspect.ratio = 1)
}

cowplot::plot_grid(plotlist = plot_list, nrow = 1, align = "v")
ggsave("beta/prok.beta_div.bray.for_region.facet_by_region.pdf", height = 3, width = 22)

cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
ggsave("beta/prok.beta_div.bray.for_region.facet_by_region.boxplot.pdf", height = 4, width = 20)

write.table(adonis2_data, "beta/prok.beta_div.bray.for_region.facet_by_region.adonis.txt", sep = "\t", quote = F, row.names = F)

adonis2_data %>% 
  mutate(name = factor(name, region_order),
         plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>% 
  ggbarplot(x = "name", y = "r2adj", label = .$plab, fill = "name", palette = colors_region,
            legend = "none", color = "black", x.text.angle = 90, lab.col = "red")
ggsave("beta/prok.beta_div.bray.for_region.facet_by_region.adonis.barplot.pdf", height = 5, width = 5)

disper_data %>% 
  write.table("beta/prok.beta_div.bray.for_region.facet_by_region.disper.txt", sep = "\t", row.names = F, quote = F)
disper_data2 %>% 
  write.table("beta/prok.beta_div.bray.for_region.facet_by_region.disper2.txt", sep = "\t", row.names = F, quote = F)

# bray prok for region facet by breed -------------------------------------

source("/code/R_func/plot_PCoA.R")
colors_region <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#F9C851","#7fc97f")
colors_breed <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f")
region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")
breed_order <- c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")

group_x <- group %>% 
  filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
  select(sample, breed = breed2, region) %>% 
  mutate(region = factor(region, region_order))

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

# bray
plot_list <- list()
plot_list2 <- list()
adonis2_data <- data.frame(name = character(), r2 = numeric(), r2adj = numeric(), pval = numeric())
disper_data <- data.frame(name = character(), pval = numeric(), method = character())
disper_data2 <- rbind()

for (i in breed_order) {
  group_i <- filter(group_x, breed == i) %>% 
    rename(group = region)
  profile_i <- select(profile, any_of(group_i$sample))
  data <- calcu_PCoA(profile_i, group = group_i, adonis2 = T, dis_method = "bray")
  
  plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line", 
                              title = paste0("Bray-distance PCoA | ", i),
                              add_lab_to_plot = T, dis_method = "bray", 
                              group_color = colors_region, show_legend = F)
  
  adonis2_data <- adonis2_data %>% 
    add_row(name = i, r2 = data$adonis2_r2, r2adj = data$adonis2_r2adj, pval = data$adonis2_p)
  
  disper <- betadisper(data$dist, data$group$group)
  
  disper_data <- disper_data %>% 
    add_row(name = i, pval = anova(disper)[1, 5], method = "anova")
  
  disper_data2 <- rbind(disper_data2, 
                        TukeyHSD(disper)[["group"]] %>% 
                          data.frame %>% 
                          rownames_to_column("group_pair") %>% 
                          add_column(group = i, method = "TukeyHSD"))
  
  plot_data <- reshape2::melt(as.matrix(data$dist)) %>% 
    filter(as.character(Var1) != as.character(Var2)) %>%
    mutate_if(is.factor, as.character) %>% 
    left_join(group_i, by = c("Var1" = "sample"))
  
  plot_list2[[i]] <- ggboxplot(plot_data, "group", "value", fill = "group", xlab = "", x.text.angle = 90,
                               ylab = "Bray curtis-based dissimilarity",
                                 legend = "none", palette = colors_region, title = paste0("Bray | ", i)) + 
    geom_line(aes(group = 1, x = group, y = median),
              plot_data %>% group_by(group) %>% summarise(median = median(value)),
              lty = "dashed") +
    stat_compare_means(aes(group = group)) +
    theme(aspect.ratio = 2/3)
}

cowplot::plot_grid(plotlist = plot_list, nrow = 1)
ggsave("beta/prok.beta_div.bray.for_region.facet_by_breed.pdf", height = 3, width = 12)

cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
ggsave("beta/prok.beta_div.bray.for_region.facet_by_breed.boxplot.pdf", height = 4, width = 15)

write.table(adonis2_data, "beta/prok.beta_div.bray.for_region.facet_by_breed.adonis.txt", 
            sep = "\t", quote = F, row.names = F)

adonis2_data %>% 
  mutate(name = factor(name, breed_order),
         plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>% 
  ggbarplot(x = "name", y = "r2adj", label = .$plab, fill = "name", palette = colors_breed,
            legend = "none", color = "black", x.text.angle = 90, lab.col = "red")
ggsave("beta/prok.beta_div.bray.for_region.facet_by_breed.adonis.barplot.pdf", height = 5, width = 5)

disper_data %>% 
  write.table("beta/prok.beta_div.bray.for_region.facet_by_breed.disper.txt", sep = "\t", row.names = F, quote = F)
disper_data2 %>% 
  write.table("beta/prok.beta_div.bray.for_region.facet_by_breed.disper2.txt", sep = "\t", row.names = F, quote = F)

# bray phage for breed -----------------------------------------------

library(patchwork)
source("/code/R_func/plot_PCoA.R")
group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/final_vOTUs.tpm.rds")
colors <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#F9C851","#ACDEF3","#f881bf",
            "#94d294","#d49264","#95beda","#eb9ecd","#F9C6B3","#fff087","#c3e6dd","#b3b3b3")

group_x <- group %>%
  filter(region == "Rectum" & group2 != "") %>% 
  select(sample, group = breed2, class) %>% 
  arrange(class, group) %>% 
  mutate(group = factor(group, unique(group)))

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

data <- calcu_PCoA(profile = profile_x, group = group_x, dis_method = "bray", adonis2 = T)

disper <- betadisper(data$dist, data$group$group) %>% 
  TukeyHSD() 
disper$group %>% 
  data.frame() %>% 
  rownames_to_column("group_pair") %>% 
  write.table("beta/phage.beta_div.bray.for_breed.disper.txt", sep = "\t", quote = F, row.names = F)

plot_data <- left_join(data$points, group_x, by = "sample")

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

ggsave("beta/phage.beta_div.bray.for_breed.manual.pdf", height = 12, width = 14)

# adonis 分析
data <- calcu_pairwise_adonis(profile_x, group_x, dis_method = "bray")

write.table(data, "beta/phage.beta_div.bray.for_breed.adonis.tsv", sep = "\t", quote = F, row.names = F)

plot_data <- data %>% 
  filter(grepl("JL-DLLW", group_pair)) %>% 
  mutate(plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", ""))),
         name = gsub("_vs_JL-DLLW", "", group_pair),
         name = factor(name, name))

ggbarplot(plot_data, x = "name", y = "r2adj", fill = "name", palette = colors, legend = "none",
          x.text.angle = 90, xlab = "", ylab = "Adjusted R2", label = plot_data$plab, 
          lab.col = "red", width = .6)

ggsave("beta/phage.beta_div.bray.for_breed.adonis.pdf", height = 6, width = 6)

levels <- levels(group_x$group)

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

ggsave("beta/phage.beta_div.bray.for_breed.adonis2.pdf", width = 6, height = 6)

# bray phage for region facet by region ------------------------------------------------------------------------------------

source("/code/R_func/plot_PCoA.R")
colors_region <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#F9C851","#7fc97f")
colors_breed <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f")
region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")

group_x <- group %>% 
  filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
  select(sample, breed = breed2, region)

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

# bray
plot_list <- list()
plot_list2 <- list()
adonis2_data <- data.frame(name = character(), r2 = numeric(), r2adj = numeric(), pval = numeric())
disper_data <- data.frame(name = character(), pval = numeric(), method = character())
disper_data2 <- rbind()

for (i in region_order) {
  group_i <- filter(group_x, region == i) %>% 
    rename(group = breed)
  profile_i <- select(profile, any_of(group_i$sample))
  data <- calcu_PCoA(profile_i, group = group_i, adonis2 = T)
  
  plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line", 
                              title = paste0("Bray-distance PCoA | ", i),
                              add_lab_to_plot = T, dis_method = "bray", 
                              group_color = colors_breed, show_legend = F)
  
  adonis2_data <- adonis2_data %>% 
    add_row(name = i, r2 = data$adonis2_r2, r2adj = data$adonis2_r2adj, pval = data$adonis2_p)
  
  disper <- betadisper(data$dist, data$group$group)
  
  disper_data <- disper_data %>% 
    add_row(name = i, pval = anova(disper)[1, 5], method = "anova")
  
  disper_data2 <- rbind(disper_data2, 
                        TukeyHSD(disper)[["group"]] %>% 
                          data.frame %>% 
                          rownames_to_column("group_pair") %>% 
                          add_column(group = i, method = "TukeyHSD"))
  
  plot_data <- reshape2::melt(as.matrix(data$dist)) %>% 
    filter(as.character(Var1) != as.character(Var2)) %>%
    mutate_if(is.factor, as.character) %>% 
    left_join(group_i, by = c("Var1" = "sample"))
  
  var <- plot_data %>% 
    filter(group == "JL-DLLW") %>% 
    pull(value) %>% 
    median()
  
  plot_list2[[i]] <- ggboxplot(plot_data, "group", "value", fill = "group", xlab = "", x.text.angle = 90,
                               ylab = "Bray curtis-based dissimilarity",
                               legend = "none", palette = colors_region, title = paste0("Bray | ", i)) + 
    stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", color = "red") +
    geom_hline(yintercept = var, lty = "dashed", linewidth = .4) +
    theme(aspect.ratio = 1)
}

cowplot::plot_grid(plotlist = plot_list, nrow = 1)
ggsave("beta/phage.beta_div.bray.for_region.facet_by_region.pdf", height = 3, width = 22)

cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
ggsave("beta/phage.beta_div.bray.for_region.facet_by_region.boxplot.pdf", height = 4, width = 20)

write.table(adonis2_data, "beta/phage.beta_div.bray.for_region.facet_by_region.adonis.txt", sep = "\t", quote = F, row.names = F)

adonis2_data %>% 
  mutate(name = factor(name, region_order),
         plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>% 
  ggbarplot(x = "name", y = "r2adj", label = .$plab, fill = "name", palette = colors_region,
            legend = "none", color = "black", x.text.angle = 90, lab.col = "red")

ggsave("beta/phage.beta_div.bray.for_region.facet_by_region.adonis.pdf", height = 5, width = 5)

disper_data %>% 
  write.table("beta/phage.beta_div.bray.for_region.facet_by_region.disper.txt", sep = "\t", row.names = F, quote = F)
disper_data2 %>% 
  write.table("beta/phage.beta_div.bray.for_region.facet_by_region.disper2.txt", sep = "\t", row.names = F, quote = F)

# bray phage for region facet by breed ------------------------------------------------------------------------------------

source("/code/R_func/plot_PCoA.R")
colors_region <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#F9C851","#7fc97f")
colors_breed <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f")
region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")
breed_order <- c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")

group_x <- group %>% 
  filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
  select(sample, breed = breed2, region) %>% 
  mutate(region = factor(region, region_order))

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

# bray
plot_list <- list()
plot_list2 <- list()
adonis2_data <- data.frame(name = character(), r2 = numeric(), r2adj = numeric(), pval = numeric())
disper_data <- data.frame(name = character(), pval = numeric(), method = character())
disper_data2 <- rbind()

for (i in breed_order) {
  group_i <- filter(group_x, breed == i) %>% 
    rename(group = region)
  profile_i <- select(profile, any_of(group_i$sample))
  data <- calcu_PCoA(profile_i, group = group_i, adonis2 = T, dis_method = "bray")
  
  plot_list[[i]] <- plot_PCoA(distance = data$dist, group = group_i, display_type = "line", 
                              title = paste0("Bray-distance PCoA | ", i),
                              add_lab_to_plot = T, dis_method = "bray", 
                              group_color = colors_region, show_legend = F)
  
  adonis2_data <- adonis2_data %>% 
    add_row(name = i, r2 = data$adonis2_r2, r2adj = data$adonis2_r2adj, pval = data$adonis2_p)
  
  disper <- betadisper(data$dist, data$group$group)
  
  disper_data <- disper_data %>% 
    add_row(name = i, pval = anova(disper)[1, 5], method = "anova")
  
  disper_data2 <- rbind(disper_data2, 
                        TukeyHSD(disper)[["group"]] %>% 
                          data.frame %>% 
                          rownames_to_column("group_pair") %>% 
                          add_column(group = i, method = "TukeyHSD"))
  
  plot_data <- reshape2::melt(as.matrix(data$dist)) %>% 
    filter(as.character(Var1) != as.character(Var2)) %>%
    mutate_if(is.factor, as.character) %>% 
    left_join(group_i, by = c("Var1" = "sample"))
  
  plot_list2[[i]] <- ggboxplot(plot_data, "group", "value", fill = "group", xlab = "", x.text.angle = 90,
                              ylab = "Bray curtis-based dissimilarity",
                              legend = "none", palette = colors_region, title = paste0("Bray | ", i)) + 
    geom_line(aes(group = 1, x = group, y = median),
              plot_data %>% group_by(group) %>% summarise(median = median(value)),
              lty = "dashed") +
    theme(aspect.ratio = 2/3)
}

cowplot::plot_grid(plotlist = plot_list, nrow = 1)
ggsave("beta/phage.beta_div.bray.for_region.facet_by_breed.pdf", height = 3, width = 12)

cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
ggsave("beta/phage.beta_div.bray.for_region.facet_by_breed.boxplot.pdf", height = 4, width = 15)

write.table(adonis2_data, "beta/phage.beta_div.bray.for_region.facet_by_breed.adonis.txt", 
            sep = "\t", quote = F, row.names = F)

adonis2_data %>% 
  mutate(name = factor(name, breed_order),
         plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>% 
  ggbarplot(x = "name", y = "r2adj", label = .$plab, fill = "name", palette = colors_breed,
            legend = "none", color = "black", x.text.angle = 90, lab.col = "red")

ggsave("beta/phage.beta_div.bray.for_region.facet_by_breed.adonis.pdf", height = 5, width = 5)

disper_data %>% 
  write.table("beta/phage.beta_div.bray.for_region.facet_by_breed.disper.txt", sep = "\t", row.names = F, quote = F)
disper_data2 %>% 
  write.table("beta/phage.beta_div.bray.for_region.facet_by_breed.disper2.txt", sep = "\t", row.names = F, quote = F)
