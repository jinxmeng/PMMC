# Jinxin Meng, 20240301, 20240703 --------------------

setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/kegg/")
pacman::p_load(dplyr, tidyr, tibble, purrr, ggplot2, ggpubr)

profile <- read.delim("rawdata/kegg.count.profile", row.names = 1)
group <- read.delim("../profile/sample_group")
KOs <- rownames(profile)

# KOs 功能注释结果 旭日图 -------------------------------------

kegg_db <- read.delim("/database/KEGG_v20230401/KO_level_A_B_C_D_Description", sep = "\t", quote = "") 

data <- kegg_db %>% 
  filter(lvD %in% KOs) %>% 
  select(-lvD, -lvDdes) %>% 
  group_by(lvA, lvAdes, lvB, lvBdes, lvC, lvCdes) %>% 
  summarise(n = n()) %>% 
  ungroup()

# 旭日图
levels <- c("lvA", "lvB", "lvC")
colors <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#F9C851","#ACDEF3","#F9C6B3")

plot_data <- data %>%
  mutate(lvA = lvAdes,
         lvB = lvBdes,
         lvC = lvCdes) %>% 
  select(lvA, lvB, lvC, n) %>% 
  group_by(lvA, lvB, lvC) %>% 
  arrange(lvA, lvB, lvC, desc(n)) %>% 
  mutate(lvA = paste0("A ", lvA),
         lvB = paste0("B ", lvB),
         lvC = paste0("C ", lvC))

plot_colors <- plot_data$lvA %>% 
  table(lvA = .) %>% 
  data.frame %>% 
  arrange(desc(Freq)) %>% 
  add_column(color = colors)

plot_colors <- plot_data %>% 
  left_join(select(plot_colors, lvA, color), by = "lvA") %>% 
  gather(key = "level", value = "name", -n, -color) %>% 
  select(name, color) %>% 
  unique() %>% 
  pull(name = name)

# write.table(plot_data, "KOs.annotations.at_various_levels.txt", sep = "\t", row.names = F, quote = F)

plot_data %>% 
  gather(key = "level", value = "name", -n) %>%
  mutate(level = factor(level, unique(level)),
         name = factor(name, unique(name))) %>% 
  arrange(level, name) %>%
  group_by(level, name) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(level) %>%
  mutate(ymax = cumsum(n),
         ymin = lag(ymax, default = 0),
         xmax = as.numeric(level),
         xmin = xmax - 1,
         xpos = (xmax + xmin)/2, 
         ypos = (ymax + ymin)/2,
         prec = ypos / 41063,
         angle = -prec * 360,
         angle = ifelse(angle < 0 & angle > -180, angle + 90, angle - 90),
         label = ifelse(n > 500, gsub("\\w__", "", name), ""),
         label = gsub("[ABC] ", "", label)) %>% 
  ggplot() +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = name),
            color = "white", linewidth = .1) +
  geom_text(aes(x = xpos, y = ypos, label = label, angle = angle), size = 2.2) +
  scale_fill_manual(values = plot_colors) +
  labs(caption = "A total of 18,087 KOs") +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "none") 
ggsave("KOs.annotations.at_various_levels.pdf", width = 9, height = 9)

# KOs 代谢功能注释结果 旭日图 -------------------------------------

levels <- c("lvB", "lvC")
colors <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#F9C851","#ACDEF3","#F9C6B3","#b3de69","#D3EDE4","#e5c494","#d9d9d9")

plot_data <- data %>%
  filter(lvAdes == "Metabolism") %>% 
  mutate(lvB = lvBdes,
         lvC = lvCdes) %>% 
  select(lvB, lvC, n) %>% 
  group_by(lvB, lvC) %>% 
  arrange(lvB, lvC, desc(n)) 

plot_colors <- plot_data$lvB %>% 
  table(lvB = .) %>% 
  data.frame %>% 
  arrange(desc(Freq)) %>% 
  add_column(color = colors)

plot_colors <- plot_data %>% 
  left_join(select(plot_colors, lvB, color), by = "lvB") %>% 
  gather(key = "level", value = "name", -n, -color) %>% 
  select(name, color) %>% 
  unique() %>% 
  pull(name = name)

write.table(plot_data, "KOs.annotations.at_various_levels.metabolism.txt", sep = "\t", row.names = F, quote = F)

plot_data %>% 
  gather(key = "level", value = "name", -n) %>%
  mutate(level = factor(level, unique(level)),
         name = factor(name, unique(name))) %>% 
  arrange(level, name) %>%
  group_by(level, name) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(level) %>%
  mutate(ymax = cumsum(n),
         ymin = lag(ymax, default = 0),
         xmax = as.numeric(level),
         xmin = xmax - 1,
         xpos = (xmax + xmin)/2, 
         ypos = (ymax + ymin)/2,
         prec = ypos / 6888,
         angle = -prec * 360,
         angle = ifelse(angle < 0 & angle > -180, angle + 90, angle - 90),
         label = ifelse(n > 50, gsub("\\w__", "", name), ""),
         label = gsub(" \\[PATH.*", "", x = label)) %>% 
  ggplot() +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = name),
            color = "white", linewidth = .1) +
  geom_text(aes(x = xpos, y = ypos, label = label, angle = angle), size = 2.2) +
  scale_fill_manual(values = plot_colors) +
  labs(caption = "A total of 4,790 KOs") +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "none") 

ggsave("KOs.annotations.at_various_levels.metabolism.pdf", width = 9, height = 9)

# KOs 数量 -------------------------------------

kegg_db <- read.delim("/database/KEGG_v20230401/KO_level_A_B_C_D_Description", sep = "\t", quote = "") 

kegg_db %>% 
  filter(lvD %in% KOs) %>% 
  select(lvAdes, lvD) %>% 
  group_by(lvAdes) %>% 
  group_modify(~unique(.x$lvD) %>% length() %>% data.frame) %>% 
  rename(value = ".") %>% 
  ungroup %>% 
  ggbarplot("lvAdes", y = "value", fill = "lvAdes", sort.val = "desc", sort.by.groups = F,
            palette = c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#F9C851","#ACDEF3","#F9C6B3"),
            x.text.angle = 30, legend = "none", label = T, xlab = "", ylab = "Number of KOs")

ggsave("KOs_annotations/KOs.annotations.lvA.count.pdf", width = 6, height = 5)

kegg_db %>% 
  filter(lvD %in% KOs & lvAdes == "Metabolism") %>%
  select(lvBdes, lvD) %>% 
  group_by(lvBdes) %>% 
  group_modify(~unique(.x$lvD) %>% length() %>% data.frame) %>% 
  rename(value = ".") %>% 
  ungroup %>% 
  filter(lvBdes != "Not included in regular maps") %>% 
  ggbarplot("lvBdes", y = "value", fill = "lvBdes", sort.val = "desc", sort.by.groups = F,
            palette = c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#F9C851",
                        "#ACDEF3","#F9C6B3","#82ccdc","#f0ad80","#b4b3d8"),
            x.text.angle = 30, legend = "none", label = T, xlab = "", ylab = "Number of KOs")

ggsave("KOs_annotations/KOs.annotations.lvB.count.pdf", width = 6, height = 5)


