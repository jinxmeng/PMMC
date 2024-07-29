# Jinxin Meng, 20240626, 20240627 ---------------------------------------

setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/diversity/")
pacman::p_load(dplyr, tidyr, tibble, purrr, ggplot2, ggpubr, ggpmisc)

library(patchwork)
source("/code/R_func/plot_PCoA.R")
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

data <- calcu_pairwise_adonis(profile_x, group_x, dis_method = "bray")

write.table(data, "beta/prok.beta_div.bray.for_breed.adonis.tsv", sep = "\t", quote = F, row.names = F)

plot_data <- data %>% 
  filter(grepl("JL-DLLW", group_pair)) %>% 
  mutate(plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", ""))),
         name = gsub("_vs_JL-DLLW", "", group_pair),
         name = factor(name, name)) %>% 
  add_row(group_pair = "", r2 = 0, r2adj = 0, pval = 0, plab = "", name = "", .before = 1) %>% 
  mutate(num = 1:nrow(.),
         angle = 90 - 360 * (num - 0.5) /nrow(.),
         hjust = ifelse( angle < -90, 1, 0),
         angle = ifelse(angle < -90, angle + 180, angle))

ggbarplot(plot_data, x = "name", y = "r2adj", fill = "name", palette = colors, legend = "none",
          x.text.angle = 90, xlab = "", ylab = "Adjusted R2", label = plot_data$plab, 
          lab.col = "red", width = .6) +
  geom_text(data = plot_data, 
            aes(x = name, y = r2adj + .1, label = name, hjust = hjust, angle = angle), 
            color="black", size = 2.5, inherit.aes = F) + 
  ylim(-0.3, 0.7) +
  coord_polar(start = 0) +
  theme(axis.line = element_blank(),
        axis.text.x = element_blank())

ggsave("prok.beta_div.bray.for_breed.adonis.circplot.pdf", height = 6, width = 6)