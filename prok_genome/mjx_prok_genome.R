# Jinxin Meng, 20240125, 20240624 -------------------------------------------------

pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr, forcats)
setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/prok_genome/")

# 分离菌基因组分析 -------------------------------------------------------------------

data <- read.delim("data.isolate.genome.metadata.txt") %>% 
  mutate(total_len_Mbp = total_len / 1e6, 
         N50_len_Kbp = log10(N50_len/1e3))

source("/code/R_func/taxa.R")
source("/code/R_func/utilities.R")

data <- read.delim("data.isolate.genome.metadata.txt") %>% 
  select(name = isolate_ID, taxonomy = gtdb_classification) %>% 
  taxa_split(rm_suffix = T)

vec <- data$family %>% table %>% sort(decreasing = T) %>% names()

plot_data <- data %>% 
  mutate(family = factor(family, vec),
         fill = forcats::fct_lump_n(family, 9, other_level = "f__Other"),
         novo = ifelse(species == "s__Unknown", "*", "")) %>% 
  arrange(fill, novo) %>% 
  add_column(x = rep(1:25, times = 13)[1:314], y = rep(13:1, each = 25)[1:314])

colors <- c("#58BBD0","#EB9256","#9B99CB","#AD9F2A","#87CCBA","#ACDEF3","#F9C6B3","#F5EAF2","#F9C851","#D3EDE4")

ggplot(plot_data, aes(x, y, fill = fill)) + 
  geom_tile(aes(height = 1, width = 1), color = "white", size = .5) + 
  geom_text(aes(label = novo)) +
  scale_fill_manual(values = colors) +
  coord_fixed(ratio = 1) +
  labs(fill = "Family level") +
  theme_void() +
  theme(legend.text = element_text(face = "italic"))

ggsave("isolate.genome_count.tile.pdf", width = 10, height = 4)

# 分离菌和MAGs聚类比较 -------------------------------------------------------------------

# 根据聚类结果模拟一个韦恩图
# genomospecies
# /share/data1/mjx/proj/08.black_pig_isolate_genomics_20231212/tmp.contrast/cls_*.res
# isolate-only 120 
# MAGs-only 3932
# both 37

library(ggvenn)

x = list("Metagenome-assembly" = 121:3992, "Isolates" = 1:157)
ggvenn(x, fill_color = c("#F9C6B3","#ACDEF3"), auto_scale = T, stroke_color = "white", stroke_size = .8, )
ggsave("isolates.cluster_contrast.genomospecies.venn.pdf", width = 5, height = 5)

# genome-based strains
# isolate-only 191
# MAGs-only 18568
# both 24

x = list("Metagenome-assembly" = 192:18568, "Isolates" = 1:215)
ggvenn(x, fill_color = c("#F9C6B3","#ACDEF3"), auto_scale = T, stroke_color = "white", stroke_size = .8)
ggsave("isolates.cluster_contrast.genomostrains.venn.pdf", width = 5, height = 5)

# 基因组菌株 统计分析 -----------------------------------------------------------------

source("/code/R_func/plot_pie.R")
data <- read.delim("data.strain.genome.metadata.txt")

data$quality2 %>% 
  table(name = .) %>% 
  data.frame %>%
  rename(n = Freq) %>% 
  plot_pie(color = "black", font_size = 4, fill = c("#82ccdc","#f0ad80","#b4b3d8"), start = 90)
ggsave("strain.genome.qs.pdf", width = 4, height = 4)

# genome size and N50
data %>% 
  mutate(x = genome_size/1e6, y = log10(N50_len/1e3)) %>% 
  ggscatterhist("x", "y", color = "source", xlab = "Genome size (Mbp)", 
              ylab = "log10 N50 (Kbp)", size = .6, margin.plot = "density", linetype = "dashed", 
              ggtheme = theme_pubr(), legend = "right")
ggsave("strain.genome_size_N50.pdf", width = 7, height = 5)

# completeness and contamination
ggscatterhist(data, x = "completeness", y = "contamination", color = "quality2", xlab = "Completeness (%)", 
              ylab = "Contamination (%)", size = .6, margin.plot = "density", linetype = "dashed", 
              ggtheme = theme_pubr(), legend = "right", palette = c("#82ccdc","#f0ad80","#b4b3d8"))
ggsave("strain.completeness_contamination.pdf", width = 9, height = 5)

# 基因组集合聚类比较 --------------------------------------------------------------------------------------------------------------------

data <- read.delim("data/cls_95.adj")  %>% 
  relocate(HuJ_2024, .after = "LiM_2024")

plot_data <- data %>%
  mutate(LiM_2024 = ifelse(LiM_2024 > 0, cluster, NA),
         HuJ_2024 = ifelse(HuJ_2024 > 0, cluster, NA),
         ChenC_2021 = ifelse(ChenC_2021 > 0, cluster, NA),
         ZhangF_2023 = ifelse(ZhangF_2023 > 0, cluster, NA)) %>% 
  column_to_rownames(var = "cluster") %>% 
  as.list() %>% 
  lapply(., function(x) as.character(na.omit(x)))

library(ggVennDiagram)
ggVennDiagram(plot_data, category.names = c("LiM_2024","HuJ_2024","ChenC_2021","ZhangF_2023"),
              label = "both", label_color = "black", label_alpha = 0, edge_lty = "dashed", edge_size = 1) +
  scale_fill_gradient(low="#fef5f5",high = "#f0ad80", name = "Number of genomospecies")

ggsave("species.contrast.venn.pdf", width = 6, height = 6)

# 2385个仅在本研究发现的基因组物种，1435个与GTDB参考基因组的ANI < 95%
data.frame(a = factor(c("unknown", "known"), c("unknown", "known")),
           b = c(1435, 2385-1435)) %>% 
  ggbarplot(y = "b", fill = "a", rotate = T, color = "white", 
            palette = c("#82ccdc","#f0ad80"), size = 1)

ggsave("species.contrast.gtdb.barplot.pdf", width = 5, height = 3)

# 稀释曲线 ------------------------------------------------------------------------

sample <- dir("rare", full.names = T)

rare_res <- list()

for (i in sample) {
  
  data <- read.delim(i, header = F, col.names = "name") %>% 
    pull()
  
  name = strsplit(i, "/")[[1]][2]
  
  rare_i <- rbind()
  for (j in c(seq(1, length(data), 200), length(data))) {
    
    rep = 1
    sampling <- list()
    while (rep <= 299) {
      sampling[[rep]] <- data[sample(1:length(data), size = j, replace = F)] %>% unique()
      rep = rep + 1
    }
    
    rare_i <- rbind(
      rare_i, 
      data.frame(name = name,
                 step = j, 
                 obs = map_vec(sampling, \(x) length(x)) %>% mean(), 
                 sd = map_vec(sampling, \(x) length(x)) %>% sd()))
  }
  
  rare_res[[name]] <- rare_i
}

reduce(rare_res, rbind) %>% 
  ggline(x = "step", y = "obs", color = "name", shape = NA, x.text.angle = 90)

ggsave("species.rare.line.pdf", width = 6, height = 5)

# 基因组物种统计 ------------------------------------------------------------------------------

source("/Code/R_func/plot_pie.R")
data <- read.delim("data.species.genome.metadata.txt") %>% 
  map_dfc(\(x) gsub("__$", "__Unknown", x))

data %>% 
  map_vec(\(x) length(unique(x[!grepl("__Unknown", x)]))) %>% 
  data.frame(n = .) %>% 
  rownames_to_column(var = "taxa") %>% 
  mutate(taxa = factor(taxa, levels = rev(taxa))) %>% 
  ggbarplot(x = "taxa", y = "n", fill = "taxa", label = T, size = .5,
            palette = c("#58BBD0","#EB9256","#9B99CB","#AD9F2A","#87CCBA","#F9C851","#ACDEF3"),
            legend = "none", xlab = "", ylab = "Number of taxa") +
  theme(aspect.ratio = 3/4)

ggsave("species.number_taxa_at_various_level.bar.pdf", width = 5, height = 4)

data %>% 
  map_vec(\(x) length(unique(x[!grepl("__Unknown", x)]))) %>% 
  data.frame(n = .) %>% 
  rownames_to_column(var = "taxa") %>% 
  ggdotchart(x = "taxa", y = "n", color = "taxa", label = T, size = 5, xtickslab.rt = 0, 
             sorting = "descending", add = "segment", ylab = "Number of taxa",
             add.params = list(color = "lightgray", size = 1), legend = "none",
             palette = c("#58BBD0","#EB9256","#9B99CB","#AD9F2A","#87CCBA","#F9C851","#ACDEF3")) +
  theme(aspect.ratio = 3/4)

ggsave("species.number_taxa_at_various_level.dot.pdf", width = 5, height = 4)

data %>% 
  map_vec(\(x) sum(grepl("Unknown", x))) %>% 
  data.frame(n = .) %>% 
  rownames_to_column(var = "taxa") %>% 
  mutate(taxa = factor(taxa, levels = rev(taxa))) %>% 
  ggbarplot(x = "taxa", y = "n", fill = "taxa", legend = "none", xlab = "", size = .5,
            ylab = "Number of genomospecies with unknown taxa", label = T,
            palette = c("#58BBD0","#EB9256","#9B99CB","#AD9F2A","#87CCBA","#F9C851","#ACDEF3")) +
  theme(aspect.ratio = 3/4)

ggsave("species.number_genomospecies_with_unknown_taxa.pdf", width = 5, height = 4)

data$phylum %>% 
  table %>% 
  data.frame %>% 
  rename(name = ".", n = "Freq") %>%
  arrange(desc(n)) %>% 
  mutate(name = factor(name, levels = name),
         name = gsub("p__", "", name)) %>% 
  ggbarplot(x = "name", y = "n", fill = "name", legend = "none", xlab = "", 
            ylab = "Number of genomospecies", label = T, x.text.angle = 90) +
  scale_fill_hue() + 
  theme(aspect.ratio = 1/3)

ggsave("species.number_genomospecies_at_phylum.bar.pdf", width = 9, height = 6)

phylum_vec <- data$phylum %>% 
  table %>% 
  data.frame(stringsAsFactors = F) %>% 
  arrange(desc(Freq)) %>% 
  pull(var = 1) %>% 
  as.character()

data %>% 
  select(name = phylum, strains) %>% 
  mutate(strains = as.numeric(strains)) %>% 
  group_by(name) %>% 
  summarise(n = sum(strains)) %>% 
  mutate(name = factor(name, levels = phylum_vec)) %>% 
  arrange(name) %>% 
  mutate(name = gsub("p__", "", name)) %>% 
  ggbarplot(x = "name", y = "n", fill = "name", legend = "none", xlab = "", 
            ylab = "Number of genome-based strains", label = T, x.text.angle = 90) +
  scale_fill_hue() + 
  theme(aspect.ratio = 1/3)

ggsave("strain.number_strain_at_phylum.bar.pdf", width = 9, height = 6)

data$family %>% 
  table %>% 
  data.frame %>% 
  rename(name = ".", n = "Freq") %>% 
  arrange(desc(n)) %>% 
  mutate(name = factor(name, levels = name)) %>% 
  plot_pie(fill = "auto", top_n = 11, font_size = 4, color = "black", start = 90)

ggsave("species.number_genomospecies_at_family.pie.pdf", width = 5, height = 5)

data$genus %>% 
  table %>%
  data.frame %>% 
  rename(name = ".", n = "Freq") %>% 
  arrange(desc(n)) %>% 
  mutate(name = factor(name, levels = name)) %>% 
  plot_pie(fill = "auto", top_n = 10, font_size = 4, color = "black", start = 90)

ggsave("species.number_genomospecies_at_genus.pie.pdf", width = 5, height = 5)

# 旭日图
taxa_levels = c("domain", "phylum", "class", "order", "family", "genus")

colors <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#F9C851","#ACDEF3","#F9C6B3","#F5EAF2","#D3EDE4")

plot_data <- data %>% 
  mutate(phylum = gsub("_[A-Z]$", "", phylum),
         genus = paste0(family, "  ", genus)) %>% 
  group_by(domain, phylum, class, order, family, genus) %>% 
  summarise(n = n()) %>% 
  arrange(domain, phylum, class, order, family, genus, desc(n)) 

plot_colors <- plot_data$phylum %>% 
  table(phylum = .) %>% 
  data.frame %>% 
  arrange(desc(Freq)) %>% 
  add_column(color = c(colors, rep("#d9d9d9", 14)))

plot_colors <- plot_data %>% 
  left_join(select(plot_colors, phylum, color), by = "phylum") %>% 
  gather(key = "level", value = "name", -n, -color) %>% 
  select(name, color) %>% 
  unique() %>% 
  pull(name = name)
plot_colors[grepl("d__Archaea", names(plot_colors))] <- "#ACDEF3"
plot_colors[grepl("d__Bacteria", names(plot_colors))] <- "#82ccdc"
# plot_colors[grepl("Unknown", names(plot_colors))] <- "#d9d9d9"

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
         prec = ypos / 4019,
         angle = -prec * 360,
         angle = ifelse(angle < 0 & angle > -180, angle + 90, angle - 90),
         label = ifelse(n > 50, gsub("\\w__", "", name) %>% gsub("\\w+  ", "", x = .), "")) %>% 
  ggplot() +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = name),
            color = "white", linewidth = .1) +
  geom_text(aes(x = xpos, y = ypos, label = label, angle = angle), size = 2.2, fontface = "italic") +
  scale_fill_manual(values = plot_colors) +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "none") 

ggsave("species.number_taxa_at_various_level.sunbrust.pdf", width = 9, height = 9)
