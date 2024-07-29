# Jinxin Meng, 20240305, 20240719 ------------------

setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/compos/")
pacman::p_load(dplyr, tidyr, tibble, purrr, ggplot2, ggpubr)
library(randomForest)

# diff function -----------------------

source("F:/Code/R_func/calcu_fisher.R")
source("F:/Code/R_func/difference_analysis.R")

meta <- read.delim("meta_analysis/meta.for_breed.genomospecies_smd.tsv") %>% 
  select(name, estimate, padj) %>% 
  distinct() %>% 
  filter(padj < 0.05) %>% 
  mutate(enriched = ifelse(estimate < 0, "CP", "NBP")) %>% 
  select(name,  enriched)

meta$enriched %>% table

data <- read.delim("diff_function/KO.tsv", header = F, col.names = c("gene", "KO")) %>% 
  mutate(name = gsub("\\.\\d+$", "", gene)) %>% 
  select(-gene) %>% 
  distinct() %>% 
  left_join(meta, by = "name") %>% 
  group_by(KO, enriched) %>% 
  summarise(value = n()) %>% 
  ungroup

data2 <- data %>% 
  spread(key = "enriched", value = "value", fill = 0) %>% 
  rename(NBP_pos = NBP, CP_pos = CP) %>% 
  mutate(NBP_neg = 2257 - NBP_pos, 
         CP_neg = 712 - CP_pos) %>% 
  select(name = KO, x_pos = "NBP_pos", y_pos = "CP_pos", x_neg = "NBP_neg", y_neg = "CP_neg")

diff <- calcu_fisher(data2, x = "NBP", y = "CP", padj = T)

# kegg pathway -----------------------

kegg_db <- read.delim("/database/KEGG_v20230401/KO_level_A_B_C_D_Description", sep = "\t", quote = "")
path_db <- kegg_db %>% 
  select(lvC, lvD) %>% 
  group_by(lvC) %>% 
  group_map(~.x %>% pull(lvD), .keep = T)
names(path_db) <- kegg_db %>% 
  group_by(lvC) %>% 
  group_keys() %>% 
  pull()

data3 <- diff %>% 
  add_plab(by = "pval") %>% 
  filter(plab != "") %>% 
  mutate(enriched = ifelse(NBP_occur > CP_occur, "NBP", "CP")) %>% 
  inner_join(kegg_db %>% select(lvD, lvAdes, lvBdes, lvC, lvCdes), by = c("name" = "lvD"))
write.table(data3, "diff_function/diff_KO.tsv", sep = "\t", row.names = F, quote = F)

data3 %>%
  # filter(lvAdes == "Metabolism") %>% 
  pull(lvAdes) %>% 
  table %>%
  data.frame() %>% 
  rename(name = ".", value = Freq) %>% 
  ggbarplot("name", "value", fill = "name", legend = "none", sort.val = "asc", sort.by.groups = F,
            xlab = "", ylab = "Number of KOs", rotate = T)
ggsave("diff_function/total_KO.lvA.stat.pdf", width = 7, height = 5)

NBP_KO <- data3 %>% 
  filter(enriched == "NBP") %>% 
  pull(1) %>% 
  unique()
CP_KO <- data3 %>% 
  filter(enriched == "CP") %>% 
  pull(1) %>% 
  unique()

NBP_ko <- map2_df(path_db, names(path_db), \(x, y) data.frame(name = y, NBP_prec = sum(NBP_KO %in% x)))
CP_ko <- map2_df(path_db, names(path_db), \(x, y) data.frame(name = y, CP_prec = sum(CP_KO %in% x)))

x <- left_join(NBP_ko, CP_ko, by = "name") %>% 
  inner_join(kegg_db %>% 
               select(lvAdes, lvBdes, lvC, lvCdes) %>% 
               distinct(), 
             by = c("name" = "lvC"))

# total
data_total <- read.delim("diff_function/KO.tsv", header = F, col.names = c("gene", "KO")) %>% 
  mutate(name = gsub("\\.\\d+$", "", gene)) %>% 
  select(-gene) %>% 
  distinct() %>% 
  left_join(meta, by = "name") %>% 
  group_by(KO, enriched) %>% 
  summarise(value = n()) %>% 
  ungroup

ko00900_KO_total <- path_db$C00900

NBP_KO <- data_total %>% 
  filter(enriched == "NBP") %>% 
  pull(KO)
CP_KO <- data_total %>% 
  filter(enriched == "CP") %>% 
  pull(KO)

m <- unique(NBP_KO, CP_KO)
m[m%in%ko00900_KO_total] %>% paste0(collapse = " ")

NBP_ko <- map2_df(path_db, names(path_db), \(x, y) data.frame(name = y, NBP_total_prec = sum(NBP_KO %in% x)))
CP_ko <- map2_df(path_db, names(path_db), \(x, y) data.frame(name = y, CP_total_prec = sum(CP_KO %in% x)))

xx <- left_join(NBP_ko, CP_ko, by = "name") %>% 
  left_join(x, by = "name") %>% 
  mutate(prec = NBP_prec/CP_prec) %>% 
  relocate(prec, .after = "CP_prec")

c(NBP_KO %in% ko00900_KO_total) %>% sum
c(CP_KO %in% ko00900_KO_total) %>% sum

part01 <- c("K00099","K01662")
part02 <- c("K12506")
part03 <- c("K06013","K13787","K24873","K00869")
part04 <- c("K01597","K11778","K01641","K12503","K00938","K00054","K00805")

p1 <- data3 %>% 
  filter(grepl("ko00900", lvCdes)) %>% 
  select(1:3, plab) %>% 
  gather(key = "class", value = "value", -name, -plab) %>% 
  mutate(class = gsub("_occur", "", class),
         value = value * 100) %>% 
  filter(name %in% part01) %>% 
  mutate(class = factor(class, c("NBP", "CP"))) %>% 
  ggbarplot("name", "value", fill = "class", palette = c("#56B4E9", "#E69F00"),
            ylab = "Occurrence (%)", position = position_dodge(width = .85))

p2 <- data3 %>% 
  filter(grepl("ko00900", lvCdes)) %>% 
  select(1:3, plab) %>% 
  gather(key = "class", value = "value", -name, -plab) %>% 
  mutate(class = gsub("_occur", "", class),
         value = value * 100) %>% 
  filter(name %in% part02) %>% 
  mutate(class = factor(class, c("NBP", "CP"))) %>% 
  ggbarplot("name", "value", fill = "class", palette = c("#56B4E9", "#E69F00"),
            ylab = "Occurrence (%)", position = position_dodge(width = .85))

p3 <- data3 %>% 
  filter(grepl("ko00900", lvCdes)) %>% 
  select(1:3, plab) %>% 
  gather(key = "class", value = "value", -name, -plab) %>% 
  mutate(class = gsub("_occur", "", class),
         value = value * 100) %>% 
  filter(name %in% part03) %>% 
  mutate(class = factor(class, c("NBP", "CP"))) %>% 
  ggbarplot("name", "value", fill = "class", palette = c("#56B4E9", "#E69F00"),
            ylab = "Occurrence (%)", position = position_dodge(width = .85))

p4 <- data3 %>% 
  filter(grepl("ko00900", lvCdes)) %>% 
  select(1:3, plab) %>% 
  gather(key = "class", value = "value", -name, -plab) %>% 
  mutate(class = gsub("_occur", "", class),
         value = value * 100) %>% 
  filter(name %in% part04) %>% 
  mutate(class = factor(class, c("NBP", "CP"))) %>% 
  ggbarplot("name", "value", fill = "class", palette = c("#56B4E9", "#E69F00"),
            ylab = "Occurrence (%)", position = position_dodge(width = .85))

library(patchwork)

(p1 + p2 + p3) / p4

ggsave("diff_function/diff_KO.pdf", width = 8, height = 7)

# host
plot_data <- read.delim("diff_function/map00900.taxa.count", header = F, col.names = c("name", "taxa", "value")) %>% 
  spread(key = "taxa", value = "value", fill = 0) %>% 
  gather(key = "taxa", value = "value", -name) %>% 
  mutate(value_log10 = log10(value),
         value_log10 = ifelse(is.infinite(value_log10), -1, value_log10))

ggscatter(plot_data, "name", "taxa", fill = "value_log10", shape = 21, color = "black", size = 4,
          x.text.angle = 90) +
  scale_fill_gradient2(low = "grey65", mid = "white", high = "#ee7e77", midpoint = 0)+
  coord_fixed()
ggsave("diff_function/map00900.taxa.count.pdf", width = 8, height = 6)
