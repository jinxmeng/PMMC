# Jinxin Meng, 20240305, 20240712 ------------------

setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/compos/")
pacman::p_load(dplyr, tidyr, tibble, purrr, ggplot2, ggpubr)

# prok NCM for breed ------------------

group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/genomospecies.tpm.b50.rds")

group_x <- group %>%
  filter(group2 != "" & region == "Rectum") %>% 
  select(sample, group = breed2, class)

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

source("F:/Code/R_func/plot_NCM.R")

breeds <- unique(group_x$group)
fit <- rbind()
plot_list <- list()
for (i in breeds) {
  data <- select(profile_x, any_of(group_x$sample[group_x$group == i]))
  data <- data[rowSums(data)!=0, colSums(data)!=0]
  fit <- rbind(fit, calcu_NCM(data, i)) 
  plot_list[[i]] <- plot_NCM(data, paste0("NCM ", i))
}
write.table(fit %>% relocate("name"), "NCM/prok.breed.NCM.fits.tsv", sep = "\t", row.names = F, quote = F)
cowplot::plot_grid(plotlist = plot_list, nrow = 3)
ggsave("NCM/prok.breed.scatter.pdf", height = 12, width = 27)

fit %>% 
  left_join(select(group_x, group, class) %>% unique, by = c("name" = "group")) %>% 
  ggplot(aes(Rsqr, mN, color = class)) +
  geom_point() +
  stat_ellipse(aes(color = class, fill = class), geom = "polygon", level = .9,
               alpha = .05, lty = 2, lwd = .3, show.legend = F) +
  scale_color_manual(values = c("#56B4E9", "#E69F00")) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
  labs(x = "Fit R square", y = "Neutral community model nM values") +
  theme_pubr() +
  theme(aspect.ratio = 1) +
  guides(color = guide_legend(position = "right"))
ggsave("NCM/prok.breed.nM_value_Rsqr.pdf", width = 5, height = 4)

# prok NCM for region ---------------------

region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")
group_x <- group %>% 
  filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
  select(sample, group, breed = breed2, region) %>% 
  mutate(region = factor(region, region_order))

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

groups <- unique(group_x$group)
fit <- rbind()
plot_list <- list()
for (i in groups) {
  data <- select(profile_x, any_of(group_x$sample[group_x$group == i]))
  data <- data[rowSums(data)!=0, colSums(data)!=0]
  try(fit <- rbind(fit, calcu_NCM(data, i)))
  try(plot_list[[i]] <- plot_NCM(data, paste0("NCM ", i)))
}
write.table(fit %>% relocate("name"), "NCM/prok.region.NCM.fits.tsv", sep = "\t", row.names = F, quote = F)
cowplot::plot_grid(plotlist = plot_list, nrow = 3)
ggsave("NCM/prok.region.scatter.pdf", height = 12, width = 28)

fit <- read.delim("NCM/prok.region.NCM.fits.tsv")

colors <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#ffed6f","#7fc97f")

fit %>% 
  left_join(select(group_x, group, breed, region) %>% unique, by = c("name" = "group")) %>% 
  ggplot(aes(Rsqr, m, color = region, shape = breed)) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf), 
            fill = "#ee7e77", alpha = .01, color = "black", linetype = "dashed") +
  geom_point() +
  scale_shape_manual(values = seq(15, 18, 1)) +
  scale_color_manual(values = colors) +
  labs(x = "Fit R square", y = "Neutral community model nM values") +
  theme_pubr() +
  theme(aspect.ratio = 1) +
  guides(color = guide_legend(position = "right"),
         shape = guide_legend(position = "right"))
ggsave("NCM/prok.region.nM_value_Rsqr.pdf", width = 5, height = 4)

# phage NCM for breed ------------------

group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/final_vOTUs.tpm.rds")

group_x <- group %>%
  filter(group2 != "" & region == "Rectum") %>% 
  select(sample, group = breed2, class)

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

breeds <- unique(group_x$group)
fit <- rbind()
plot_list <- list()
for (i in breeds) {
  data <- select(profile_x, any_of(group_x$sample[group_x$group == i]))
  data <- data[rowSums(data)!=0, colSums(data)!=0]
  fit <- rbind(fit, calcu_NCM(data, i)) 
  plot_list[[i]] <- plot_NCM(data, paste0("NCM ", i))
}
write.table(fit %>% relocate("name"), "NCM/phage.breed.NCM.fits.tsv", sep = "\t", row.names = F, quote = F)
cowplot::plot_grid(plotlist = plot_list, nrow = 3)
ggsave("NCM/phage.breed.scatter.pdf", height = 12, width = 27)

fit <- read.delim("NCM/phage.breed.NCM.fits.tsv")

fit %>% 
  left_join(select(group_x, group, class) %>% unique, by = c("name" = "group")) %>% 
  ggplot(aes(Rsqr, mN, color = class)) +
  geom_point() +
  stat_ellipse(aes(color = class, fill = class), geom = "polygon", level = .9,
               alpha = .05, lty = 2, lwd = .3, show.legend = F) +
  scale_color_manual(values = c("#56B4E9", "#E69F00")) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
  labs(x = "Fit R square", y = "Neutral community model nM values") +
  theme_pubr() +
  theme(aspect.ratio = 1) +
  guides(color = guide_legend(position = "right"))
ggsave("NCM/phage.breed.nM_value_Rsqr.pdf", width = 5, height = 4)

# phage NCM for region ---------------------

region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")
group_x <- group %>% 
  filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
  select(sample, group, breed = breed2, region) %>% 
  mutate(region = factor(region, region_order))

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

groups <- unique(group_x$group)
fit <- rbind()
plot_list <- list()
for (i in groups) {
  data <- select(profile_x, any_of(group_x$sample[group_x$group == i]))
  data <- data[rowSums(data)!=0, colSums(data)!=0]
  try(fit <- rbind(fit, calcu_NCM(data, i)))
  try(plot_list[[i]] <- plot_NCM(data, paste0("NCM ", i)))
}
write.table(fit %>% relocate("name"), "NCM/phage.region.NCM.fits.tsv", sep = "\t", row.names = F, quote = F)
cowplot::plot_grid(plotlist = plot_list, nrow = 3)
ggsave("NCM/phage.region.scatter.pdf", height = 12, width = 28)

fit <- read.delim("NCM/phage.region.NCM.fits.tsv")

colors <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#ffed6f","#7fc97f")

fit %>% 
  left_join(select(group_x, group, breed, region) %>% unique, by = c("name" = "group")) %>% 
  ggplot(aes(Rsqr, mN, color = region, shape = breed)) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf), 
            fill = "#ee7e77", alpha = .01, color = "black", linetype = "dashed") +
  geom_point() +
  scale_shape_manual(values = seq(15, 18, 1)) +
  scale_color_manual(values = colors) +
  labs(x = "Fit R square", y = "Neutral community model nM values") +
  theme_pubr() +
  theme(aspect.ratio = 1) +
  guides(color = guide_legend(position = "right"),
         shape = guide_legend(position = "right"))
ggsave("NCM/phage.region.nM_value_Rsqr.pdf", width = 5, height = 4)
