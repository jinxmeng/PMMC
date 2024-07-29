# Jinxin Meng, 20240301, 20240709 ---------------------------------------

setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/diversity/")
pacman::p_load(dplyr, tidyr, tibble, purrr, ggplot2, ggpubr, ggpmisc)

# prok for breed -----------------------------------------------

source("/code/R_func/diversity.R")
source("/code/R_func/calcu_metafor.R")
group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/genomospecies.tpm.b50.rds")

group_x <- group %>%
  filter(region == "Rectum" & group2 != "") %>% 
  select(sample, group = breed2, class) %>% 
  arrange(class, group) %>% 
  mutate(group = factor(group, unique(group)))

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

# shannon
data <- calcu_alpha(profile_x, method = "shannon") %>% 
  left_join(group_x, by = "sample")

var <- data %>% 
  filter(grepl("PN5R", sample)) %>% 
  pull(2) %>% 
  median()

p1 <- ggboxplot(data, x = "group", y = "val", xlab = "", ylab = "Shannon Index", 
                fill = "class", legend = "none", rotate = T, color = "black",
                outlier.shape = NA, palette = c("#56B4E9", "#E69F00")) +
  stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", hide.ns = T, color = "red") +
  geom_hline(yintercept = var, lty = "dashed", linewidth = .4) +
  theme(aspect.ratio = 2)

vec <- data %>% 
  filter(group == "JL-DLLW") %>% 
  pull(val)

data1 <- data %>% 
  select(-class, -sample) %>% 
  filter(group != "JL-DLLW") %>% 
  group_by(group) %>% 
  summarise(d_mean = mean(val), d_sd = sd(val), d_n = n()) %>% 
  ungroup %>% 
  add_column(c_mean = mean(vec),
             c_sd = sd(vec),
             c_n = length(vec)) %>% 
  rename(name = group) %>% 
  metafor_fit.1() %>% 
  add_column(index = "Shannon")

# richness
data <- calcu_alpha(profile_x, method = "richness") %>% 
  left_join(group_x, by = "sample")

var <- data %>% 
  filter(grepl("PN5R", sample)) %>% 
  pull(2) %>% 
  median()

p2 <- ggboxplot(data, x = "group", y = "val", xlab = "", ylab = "Richness Index", 
                fill = "class", legend = "none", rotate = T, color = "black",
                outlier.shape = NA, palette = c("#56B4E9", "#E69F00")) +
  stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", hide.ns = T, color = "red") +
  geom_hline(yintercept = var, lty = "dashed", linewidth = .4) +
  theme(aspect.ratio = 2)

vec <- data %>% 
  filter(group == "JL-DLLW") %>% 
  pull(val)

data2 <- data %>% 
  select(-class, -sample) %>% 
  filter(group != "JL-DLLW") %>% 
  group_by(group) %>% 
  summarise(d_mean = mean(val), d_sd = sd(val), d_n = n()) %>% 
  ungroup %>% 
  add_column(c_mean = mean(vec),
             c_sd = sd(vec),
             c_n = length(vec)) %>% 
  rename(name = group) %>% 
  metafor_fit.1() %>% 
  add_column(index = "Richness")

# pd
tr <- ape::read.tree("../prok_genome/tree/mjx_black_pigs_4019_species_20240706.tree.nwk")

data <- calcu_alpha(profile_x, method = "pd", tree = tr) %>% 
  left_join(group_x, by = "sample")

var <- data %>% 
  filter(grepl("PN5R", sample)) %>% 
  pull(2) %>% 
  median()

p3 <- ggboxplot(data, x = "group", y = "val", xlab = "", ylab = "Faith's Phylogenetic Diversity",
                fill = "class", legend = "none", rotate = T, color = "black",
                outlier.shape = NA, palette = c("#56B4E9", "#E69F00")) +
  stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", hide.ns = T, color = "red") +
  geom_hline(yintercept = var, lty = "dashed", linewidth = .4) +
  theme(aspect.ratio = 2)

vec <- data %>% 
  filter(group == "JL-DLLW") %>% 
  pull(val)

data3 <- data %>% 
  select(-class, -sample) %>% 
  filter(group != "JL-DLLW") %>% 
  group_by(group) %>% 
  summarise(d_mean = mean(val), d_sd = sd(val), d_n = n()) %>% 
  ungroup %>% 
  add_column(c_mean = mean(vec),
             c_sd = sd(vec),
             c_n = length(vec)) %>% 
  rename(name = group) %>% 
  metafor_fit.1() %>% 
  add_column(index = "Faith's Phylogenetic Diversity")

cowplot::plot_grid(p1, p2, p3, nrow = 1)
ggsave("alpha/prok.alpha_div.for_breed.pdf", height = 5, width = 12)

# meta-analysis
plot_data <- rbind(data1, data2, data3) %>% 
  mutate(plab = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", ""))))

write.table(plot_data, "alpha/prok.alpha_div.for_breed.meta-analysis.tsv", sep = "\t", quote = F, row.names = F)

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

ggsave("alpha/prok.alpha_div.for_breed.meta-analysis.pdf", height = 2, width = 6)

# prok for region -------------------------------

source("F:/code/R_func/diversity.R")
source("/code/R_func/calcu_metafor.R")
group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/genomospecies.tpm.b50.rds")

group_x <- group %>% 
  filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
  select(sample, breed = breed2, region)

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")

# shannon
data <- calcu_alpha(profile_x, method = "shannon") %>% 
  merge(group_x, by = "sample") %>% 
  mutate(region = factor(region, region_order))

plot_data <- data %>% 
  group_by(breed, region) %>% 
  summarise(mean = mean(val), sd = sd(val), n = n()) %>% 
  mutate(se = sd / sqrt(n))

p1 <- ggerrorplot(plot_data, x = "region", y = "mean", color = "breed", 
                  position = position_dodge(width = 0), xlab = "Location", 
                  ylab = "Shannon Index", x.text.angle = 90) +
  geom_ribbon(aes(x = region, ymin = mean - se, ymax = mean + se, fill = breed, group = breed),
              plot_data, show.legend = F, inherit.aes = F, alpha = .3) +
  stat_compare_means(aes(x = region, y = val, group = breed, color = breed), 
                     data, inherit.aes = F,
                     label = "p.signif", color = "red") +
  scale_color_manual(values = c("#82ccdc","#b4b3d8","#c2b75f","#F9C851")) +
  scale_fill_manual(values = c("#82ccdc","#b4b3d8","#c2b75f","#F9C851")) +
  theme(aspect.ratio = 1)

# richness
data <- calcu_alpha(profile_x, method = "richness") %>% 
  merge(group_x, by = "sample") %>% 
  mutate(region = factor(region, region_order))

plot_data <- data %>% 
  group_by(breed, region) %>% 
  summarise(mean = mean(val), sd = sd(val), n = n()) %>% 
  mutate(se = sd / sqrt(n))

p2 <- ggerrorplot(plot_data, x = "region", y = "mean", color = "breed", 
                  position = position_dodge(width = 0), xlab = "Location", 
                  ylab = "Richness Index", x.text.angle = 90) +
  geom_ribbon(aes(x = region, ymin = mean - se, ymax = mean + se, fill = breed, group = breed),
              plot_data, show.legend = F, inherit.aes = F, alpha = .3) +
  stat_compare_means(aes(x = region, y = val, group = breed, color = breed), 
                     data, inherit.aes = F,
                     label = "p.signif", color = "red") +
  scale_color_manual(values = c("#82ccdc","#b4b3d8","#c2b75f","#F9C851")) +
  scale_fill_manual(values = c("#82ccdc","#b4b3d8","#c2b75f","#F9C851")) +
  theme(aspect.ratio = 1)

# pd
tr <- ape::read.tree("../prok_genome/tree/mjx_black_pigs_4019_species_20240706.tree.nwk")

data <- calcu_alpha(profile_x, method = "pd", tree = tr) %>% 
  merge(group_x, by = "sample") %>% 
  mutate(region = factor(region, region_order))

plot_data <- data %>% 
  group_by(breed, region) %>% 
  summarise(mean = mean(val), sd = sd(val), n = n()) %>% 
  mutate(se = sd / sqrt(n))

p3 <- ggerrorplot(plot_data, x = "region", y = "mean", color = "breed", 
                  position = position_dodge(width = 0), xlab = "Location", 
                  ylab = "Faith's Phylogenetic Diversity Index", x.text.angle = 90) +
  geom_ribbon(aes(x = region, ymin = mean - se, ymax = mean + se, fill = breed, group = breed),
              plot_data, show.legend = F, inherit.aes = F, alpha = .3) +
  stat_compare_means(aes(x = region, y = val, group = breed, color = breed), 
                     data, inherit.aes = F, 
                     label = "p.signif", color = "red") +
  scale_color_manual(values = c("#82ccdc","#b4b3d8","#c2b75f","#F9C851")) +
  scale_fill_manual(values = c("#82ccdc","#b4b3d8","#c2b75f","#F9C851")) +
  theme(aspect.ratio = 1)

cowplot::plot_grid(p1, p2, p3, nrow = 1, align = "v")
ggsave("alpha/prok.alpha_div.for_region.pdf", height = 6, width = 12)

# phage for breed -----------------------------------------------

source("/code/R_func/diversity.R")
source("/code/R_func/calcu_metafor.R")
group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/final_vOTUs.tpm.rds")

group_x <- group %>%
  filter(region == "Rectum" & group2 != "") %>% 
  select(sample, group = breed2, class) %>% 
  arrange(class, group) %>% 
  mutate(group = factor(group, unique(group)))

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

# shannon
data <- calcu_alpha(profile_x, method = "shannon") %>% 
  left_join(group_x, by = "sample")

var <- data %>% 
  filter(grepl("PN5R", sample)) %>% 
  pull(2) %>% 
  median()

p1 <- ggboxplot(data, x = "group", y = "val", xlab = "", ylab = "Shannon Index", 
                fill = "class", legend = "none", rotate = T, color = "black",
                outlier.shape = NA, palette = c("#56B4E9", "#E69F00")) +
  stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", hide.ns = T, color = "red") +
  geom_hline(yintercept = var, lty = "dashed", linewidth = .4) +
  theme(aspect.ratio = 2)

vec <- data %>% 
  filter(group == "JL-DLLW") %>% 
  pull(val)

data1 <- data %>% 
  select(-class, -sample) %>% 
  filter(group != "JL-DLLW") %>% 
  group_by(group) %>% 
  summarise(d_mean = mean(val), d_sd = sd(val), d_n = n()) %>% 
  ungroup %>% 
  add_column(c_mean = mean(vec),
             c_sd = sd(vec),
             c_n = length(vec)) %>% 
  rename(name = group) %>% 
  metafor_fit.1() %>% 
  add_column(index = "Shannon")

# richness
data <- calcu_alpha(profile_x, method = "richness") %>% 
  left_join(group_x, by = "sample")

var <- data %>% 
  filter(grepl("PN5R", sample)) %>% 
  pull(2) %>% 
  median()

p2 <- ggboxplot(data, x = "group", y = "val", xlab = "", ylab = "Richness Index", 
                fill = "class", legend = "none", rotate = T, color = "black",
                outlier.shape = NA, palette = c("#56B4E9", "#E69F00")) +
  stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", hide.ns = T, color = "red") +
  geom_hline(yintercept = var, lty = "dashed", linewidth = .4) +
  theme(aspect.ratio = 2)

vec <- data %>% 
  filter(group == "JL-DLLW") %>% 
  pull(val)

data2 <- data %>% 
  select(-class, -sample) %>% 
  filter(group != "JL-DLLW") %>% 
  group_by(group) %>% 
  summarise(d_mean = mean(val), d_sd = sd(val), d_n = n()) %>% 
  ungroup %>% 
  add_column(c_mean = mean(vec),
             c_sd = sd(vec),
             c_n = length(vec)) %>% 
  rename(name = group) %>% 
  metafor_fit.1() %>% 
  add_column(index = "Richness")

cowplot::plot_grid(p1, p2, nrow = 1)
ggsave("alpha/phage.alpha_div.for_breed.pdf", height = 5, width = 9)

# meta-analysis
plot_data <- rbind(data1, data2) %>% 
  mutate(plab = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", ""))))

write.table(plot_data, "alpha/phage.alpha_div.for_breed.meta-analysis.tsv", sep = "\t", quote = F, row.names = F)

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
            plot_data, color = "red") +
  labs(y = "", x = "Standradized Mean Difference (Random Effect Model)") +
  scale_color_viridis_d(begin = .4) +
  theme_classic2() +
  theme(axis.text = element_text(color = "black", size = 8),
        axis.title = element_text(size = 8),
        axis.line = element_line(linewidth = .4, color = "#000000"),
        axis.ticks = element_line(linewidth = .4, color = "#000000"))

ggsave("alpha/phage.alpha_div.for_breed.meta-analysis.pdf", height = 1.6, width = 6)

# phage for region -------------------------------------------

source("/code/R_func/diversity.R")
source("/code/R_func/calcu_metafor.R")
group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/final_vOTUs.tpm.rds")

group_x <- group %>% 
  filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
  select(sample, breed = breed2, region)

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")

# shannon
data <- calcu_alpha(profile_x, method = "shannon") %>% 
  merge(group_x, by = "sample") %>% 
  mutate(region = factor(region, region_order))

plot_data <- data %>% 
  group_by(breed, region) %>% 
  summarise(mean = mean(val), sd = sd(val), n = n()) %>% 
  mutate(se = sd / sqrt(n))

p1 <- ggerrorplot(plot_data, x = "region", y = "mean", color = "breed", 
                  position = position_dodge(width = 0), xlab = "Location", 
                  ylab = "Shannon Index", x.text.angle = 90) +
  geom_ribbon(aes(x = region, ymin = mean - se, ymax = mean + se, fill = breed, group = breed),
              plot_data, show.legend = F, inherit.aes = F, alpha = .3) +
  stat_compare_means(aes(x = region, y = val, group = breed, color = breed), 
                     data, inherit.aes = F,
                     label = "p.signif", color = "red") +
  scale_color_manual(values = c("#82ccdc","#b4b3d8","#c2b75f","#F9C851")) +
  scale_fill_manual(values = c("#82ccdc","#b4b3d8","#c2b75f","#F9C851")) +
  theme(aspect.ratio = 1)

# richness
data <- calcu_alpha(profile_x, method = "richness") %>% 
  merge(group_x, by = "sample") %>% 
  mutate(region = factor(region, region_order))

plot_data <- data %>% 
  group_by(breed, region) %>% 
  summarise(mean = mean(val), sd = sd(val), n = n()) %>% 
  mutate(se = sd / sqrt(n))

p2 <- ggerrorplot(plot_data, x = "region", y = "mean", color = "breed", 
                  position = position_dodge(width = 0), xlab = "Location", 
                  ylab = "Richness Index", x.text.angle = 90) +
  geom_ribbon(aes(x = region, ymin = mean - se, ymax = mean + se, fill = breed, group = breed),
              plot_data, show.legend = F, inherit.aes = F, alpha = .3) +
  stat_compare_means(aes(x = region, y = val, group = breed, color = breed), 
                     data, inherit.aes = F,
                     label = "p.signif", color = "red") +
  scale_color_manual(values = c("#82ccdc","#b4b3d8","#c2b75f","#F9C851")) +
  scale_fill_manual(values = c("#82ccdc","#b4b3d8","#c2b75f","#F9C851")) +
  theme(aspect.ratio = 1)

cowplot::plot_grid(p1, p2, nrow = 1, align = "v")
ggsave("alpha/phage.alpha_div.for_region.pdf", height = 6, width = 8)
