# Jinxin Meng, 20240521, 20240716 --------------------

setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/BGCs/")
pacman::p_load(dplyr, tidyr, tibble, purrr, ggplot2, ggpubr, ggpmisc)

# BGCs class pie ---------------------

source("/code/R_func/plot_pie.R")
data <- read.delim("BGCs.metadata.tsv") %>% 
  select(name = BiGSCAPE_class) %>% 
  group_by(name) %>% 
  summarise(n = n())
plot_pie(data, fill = "auto", add_n = T, font_size = 4, start = 60, title = "15,919 BGCs",
         color = "black")
ggsave("BGCs.class.pie.pdf", width = 4, height = 4)

# BGCs class & phylum circular barplot -------------------------------------------

color_phylum <- structure(c("#bebada","#ffffb3","#8dd3c7","#FF9D9A","#fee08b","#a7c9e0","#de77ae","#fabfd2","#b8e186","#d9d9d9"),
                          names = c("Actinomycetota","Bacillota","Bacteroidota","Cyanobacteriota","Methanobacteriota", "Other","Pseudomonadota","Spirochaetota","Thermoproteota","Verrucomicrobiota"))

taxa <- read.delim("../prok_genome/data.strain.genome.metadata.txt")
BGCs <- read.delim("rawdata/Network_Annotations_Full.tsv", check.names = F)

data <- BGCs %>% 
  select(BGC, class = "BiG-SCAPE class") %>% 
  filter(!grepl("BGC", BGC)) %>% 
  mutate(genome = ifelse(grepl("BGC", x = BGC), "", stringr::str_extract(BGC, "[A-Za-z0-9.]+"))) %>% 
  left_join(taxa %>% select(name, phylum), by = c("genome" = "name")) %>% 
  mutate(phylum = gsub("p__", "", phylum),
         phylum = gsub("_\\w$", "", phylum),
         phylum = ifelse(phylum %in% names(color_phylum), phylum, "Other")) %>% 
  select(phylum, class) %>% 
  group_by(phylum, class) %>% 
  summarise(n = n()) %>% 
  mutate(n_log10 = log10(n), 
         class = factor(class),
         n_lab = prettyNum(n, big.mark = ","))

# plot_dat
empty_bar <- 2
empty_df <- data.frame(matrix(NA, empty_bar * nlevels(data$class), ncol(data))) %>% 
  rename(any_of(structure(colnames(.), names = colnames(data))))
empty_df$class <- rep(levels(data$class), each = empty_bar)
plot_data <- rbind(data, empty_df) %>% 
  arrange(class, phylum) %>% 
  add_column(id = 1:nrow(.))

# plot_lab
plot_label <- plot_data %>% 
  mutate(cols = nrow(.),
         angle = 90 - 360 * (id - 0.5) / cols,
         hjust = ifelse(angle < -90, 1, 0),
         angle = ifelse(angle < -90, angle + 180, angle))

plot_title <- plot_data %>% 
  group_by(class) %>% 
  summarize(min = min(id), max = max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(mean = mean(c(min, max)))

plot_grid <- plot_data %>% 
  filter(is.na(n)) %>% 
  select(class, id) %>% 
  group_by(class) %>% 
  summarise(min = min(id), max = max(id))

ggplot(plot_data, aes(x = as.factor(id), y = n_log10, fill = class)) +
  geom_bar(stat = "identity", alpha = 0.95) +
  ylim(-3, max(na.omit(plot_data$n_log10)) + 1) +
  coord_polar() +
  theme_void() +
  geom_segment(data = plot_grid, aes(x = min, y = 1, xend = max, yend = 1), color = "grey", alpha = 1, size = 0.3, inherit.aes = F) +
  geom_segment(data = plot_grid, aes(x = min, y = 2, xend = max, yend = 2), color = "grey", alpha = 1, size = 0.3, inherit.aes = F) +
  geom_segment(data = plot_grid, aes(x = min, y = 3, xend = max, yend = 3), color = "grey", alpha = 1, size = 0.3, inherit.aes = F) +
  annotate("text", x = rep(.5, 3), y = c(1, 2, 3), label = c("1", "2", "3"), color = "grey", size = 3, angle = 0, fontface = "bold", hjust = 1) +
  geom_text(data = plot_label, aes(x = id, y = n_log10 + .2, label = phylum, hjust = hjust, angle = angle), color = "black", fontface = "plain", size = 2.5, inherit.aes = F) +
  geom_text(data = plot_label, aes(x = id, y = n_log10 - .2, label = n_lab, hjust = -hjust, angle = angle), color = "black", fontface = "plain", size = 2.5, inherit.aes = F) +
  geom_segment(data = plot_title, aes(x = min, y = -.3, xend = max, yend = -.3), color = "black", size = 0.8, inherit.aes = F) +
  geom_text(data = plot_title, aes(x = mean, y = -.5, label = class, color = class), hjust = .5, vjust = .5, alpha = 0.95, size = 4, fontface = "bold", inherit.aes = F) +
  scale_fill_manual(values = c("#82ccdc","#f0ad80","#b4b3d8","#77b0c5","#c2b75f","#87CCBA","#FEDFB2","#f08b85")) +
  scale_color_manual(values = c("#82ccdc","#f0ad80","#b4b3d8","#77b0c5","#c2b75f","#87CCBA","#FEDFB2","#f08b85")) +
  guides(fill = "none", color = "none")
ggsave("BGCs.class.cirbar.pdf", width = 7, height = 7)

# BGCs class & phylum bar -------------------

plot_data <- BGCs %>% 
  select(BGC, class = "BiG-SCAPE class") %>% 
  filter(!grepl("BGC", BGC)) %>% 
  mutate(genome = ifelse(grepl("BGC", x = BGC), "", stringr::str_extract(BGC, "[A-Za-z0-9.]+"))) %>% 
  left_join(taxa %>% select(name, phylum), by = c("genome" = "name")) %>% 
  mutate(phylum = gsub("p__", "", phylum),
         phylum = gsub("_\\w$", "", phylum),
         phylum = ifelse(phylum %in% names(color_phylum), phylum, "Other")) %>% 
  group_by(phylum, class) %>% 
  summarise(n = n())

phylum_order <- aggregate(n~phylum, plot_data, sum) %>%
  arrange(desc(n)) %>% 
  select(phylum) %>% 
  pull()

plot_data <- plot_data %>% 
  mutate(phylum = factor(phylum, phylum_order)) %>% 
  group_by(phylum) %>% 
  group_modify(~.x %>% mutate(perc = n / sum(n) * 100)) %>% 
  ungroup

write.table(plot_data, "BGCs.class.phylum.bar.tsv", sep = "\t", row.names = F, quote = F)

plot_label <- aggregate(n~phylum, plot_data, sum) %>% 
  mutate(lab = prettyNum(n, big.mark = ","))
  
p1 <- ggbarplot(plot_data, "phylum", "perc", fill = "class", rotate = T, ylab = "BGCs frequency",
                palette = c("#82ccdc","#f0ad80","#b4b3d8","#77b0c5","#c2b75f","#87CCBA","#FEDFB2","#f08b85")) +
  theme(aspect.ratio = 1.5)

p2 <- ggbarplot(plot_data, "phylum", "n", fill = "class", rotate = T, xlab = "", ylab = "Number of BGCs",
                palette = c("#82ccdc","#f0ad80","#b4b3d8","#77b0c5","#c2b75f","#87CCBA","#FEDFB2","#f08b85")) +
  geom_text(data = plot_label, aes(x = phylum, y = n + 12, label = lab), hjust = 0) +
  theme(aspect.ratio = 1.5,
        axis.text.y = element_blank())

cowplot::plot_grid(p1, p2, align = "v")
ggsave("BGCs.class.phylum.bar.pdf", width = 10, height = 8)

# BGCs class & family bar -------------------

plot_data <- BGCs %>% 
  select(BGC, class = "BiG-SCAPE class") %>% 
  filter(!grepl("BGC", BGC)) %>% 
  mutate(genome = ifelse(grepl("BGC", x = BGC), "", stringr::str_extract(BGC, "[A-Za-z0-9.]+"))) %>% 
  left_join(taxa %>% select(name, family), by = c("genome" = "name")) %>% 
  mutate(family = gsub("f__", "", family),
         family = gsub("_\\w$", "", family),
         family = forcats::fct_lump_n(family, 29, other_level = "Other")) %>% 
  group_by(family, class) %>% 
  summarise(n = n())

family_order <- aggregate(n~family, plot_data, sum) %>%
  arrange(desc(n)) %>% 
  select(family) %>% 
  pull()

plot_data <- plot_data %>% 
  mutate(family = factor(family, family_order)) %>% 
  group_by(family) %>% 
  group_modify(~.x %>% mutate(perc = n / sum(n) * 100)) %>% 
  ungroup

write.table(plot_data, "BGCs.class.family.bar.tsv", sep = "\t", row.names = F, quote = F)

plot_label <- aggregate(n~family, plot_data, sum) %>% 
  mutate(lab = prettyNum(n, big.mark = ","))

p1 <- ggbarplot(plot_data, "family", "perc", fill = "class", rotate = T, ylab = "BGCs frequency",
                palette = c("#82ccdc","#f0ad80","#b4b3d8","#77b0c5","#c2b75f","#87CCBA","#FEDFB2","#f08b85")) +
  theme(aspect.ratio = 2.5)

p2 <- ggbarplot(plot_data, "family", "n", fill = "class", rotate = T, xlab = "", ylab = "Number of BGCs",
                palette = c("#82ccdc","#f0ad80","#b4b3d8","#77b0c5","#c2b75f","#87CCBA","#FEDFB2","#f08b85")) +
  geom_text(data = plot_label, aes(x = family, y = n + 12, label = lab), hjust = 0) +
  theme(aspect.ratio = 2.5,
        axis.text.y = element_blank())

cowplot::plot_grid(p1, p2, align = "v")
ggsave("BGCs.class.family.bar.pdf", width = 15, height = 15)

# BGCs class stat --------------------------------------------------------------

plot_data <- BGCs %>% 
  select(BGC, class = "BiG-SCAPE class") %>% 
  filter(!grepl("BGC", BGC)) %>% 
  mutate(genome = ifelse(grepl("BGC", x = BGC), "", stringr::str_extract(BGC, "[A-Za-z0-9.]+"))) %>% 
  left_join(taxa %>% select(name, phylum), by = c("genome" = "name")) %>% 
  mutate(phylum = gsub("p__", "", phylum),
         phylum = gsub("_\\w$", "", phylum),
         phylum = ifelse(phylum %in% names(color_phylum), phylum, "Other")) %>% 
  group_by(phylum, class) %>% 
  summarise(n = n())

plot_data <- BGCs %>% 
  select(BGC, class = "BiG-SCAPE class") %>% 
  filter(!grepl("BGC", BGC)) %>% 
  mutate(genome = ifelse(grepl("BGC", x = BGC), "", stringr::str_extract(BGC, "[A-Za-z0-9.]+"))) %>% 
  left_join(taxa %>% select(name, phylum, family, genus), by = c("genome" = "name")) %>% 
  mutate(phylum = gsub("_\\w$", "", phylum),
         family = ifelse(family == "f__", "f__Unknown", family),
         family = gsub("_\\w$", "", family),
         genus = ifelse(genus == "g__", "g__Unknown", genus),
         genus = gsub("_\\w$", "", genus)) %>% 
  group_by(phylum, family, genus, class) %>% 
  summarise(n = n())
write.table(plot_data, "BGCs.class.taxa.stat.tsv", sep = "\t", row.names = F, quote = F)

plot_data <- BGCs %>% 
  select(BGC, class = "BiG-SCAPE class") %>% 
  filter(!grepl("BGC", BGC)) %>% 
  mutate(genome = ifelse(grepl("BGC", x = BGC), "", stringr::str_extract(BGC, "[A-Za-z0-9.]+"))) %>% 
  left_join(taxa %>% select(name, phylum, family), by = c("genome" = "name")) %>% 
  mutate(phylum = gsub("_\\w$", "", phylum),
         family = ifelse(family == "f__", "f__Unknown", family),
         family = gsub("_\\w$", "", family)) %>% 
  group_by(phylum, family, class) %>% 
  summarise(n = n())
write.table(plot_data, "BGCs.class.taxa.stat2.tsv", sep = "\t", row.names = F, quote = F)

# BGCs production other --------------------------------------------------------------

metadata <- read.delim("BGCs.metadata.tsv")

data <- strsplit(metadata$Product_Prediction, "\\.") %>% 
  map_df(\(x) data.frame(name = x, n = 1/length(x))) %>% 
  mutate(name = sub("-class.*", "", name)) %>% 
  group_by(name) %>% 
  summarise(n = sum(n)) %>% 
  mutate(n = round(n, 0))
write.table(data, "BGCs.product.pieplot.txt", row.names = F, quote = F)
plot_pie(data, top_n = 13, fill = "auto", add_n = T, start = 110, font_size = 4, color = "black")
ggsave("BGCs.product.pieplot.pdf", width = 4, height = 4)

# BGCs MIBiG GCF intersect --------------------------------------------------------------------

metadata <- read.delim("BGCs.metadata.tsv")

files <- list.files("bigscape.out/network_files/2024-05-21_12-23-16_hybrids_glocal/", 
                    full.names = T, recursive = T, pattern = "0.30_0.70.tsv")

data <- map_dfr(files, \(x) read.delim(x, check.names = F)) %>% 
  select(BGC = 1, GCC = 2, GCF = 3) %>% 
  mutate(GCC = paste0("GCC", sprintf("%05d", GCC)),
         GCF = paste0("GCF", sprintf("%05d", GCF)),
         class = ifelse(grepl("BGC", BGC), "MIBiG", "LiM_2024"))

BGC_data <- data %>% 
  filter(grepl("BGC", x = BGC))

data$GCF %>% 
  unique() %>% 
  length

plot_data <- data %>% 
  group_by(GCF) %>% 
  group_modify(~.x %>% 
                 group_by(class) %>% 
                 summarise(n = n())) 

plot_data %>% 
  spread(key = "class", value = "n", fill = 0) %>% 
  column_to_rownames(var = "GCF") %>%
  apply(2, \(x) ifelse(x > 0, rownames(.), NA)) %>% 
  data.frame %>% 
  map(\(x) as.character(na.omit(x))) %>% 
  ggvenn::ggvenn(., fill_color = c("#f1b389", "#8ad0de"), stroke_linetype = "dashed", 
                 auto_scale = T, fill_alpha = .9)
ggsave("GCFs.vs.MIBiG.venn.pdf", width = 5, height = 4)

# net intersect
plot_data <- data %>% 
  group_by(GCF) %>% 
  group_modify(~.x %>% 
                 group_by(class) %>% 
                 summarise(n = n())) 

intersect <- plot_data %>% 
  spread(key = "class", value = "n", fill = 0) %>% 
  column_to_rownames(var = "GCF") %>%
  filter(MIBiG > 0 & LiM_2024 > 0) %>% 
  rownames()

files <- list.files("bigscape.out/network_files/2024-05-21_12-23-16_hybrids_glocal/", 
                    full.names = T, recursive = T, pattern = "30.network")

network <- map_dfr(files, \(x) read.delim(x, check.names = F)) %>% 
  select(from = 1, to = 2, dist = 3)

data <- data %>% 
  select(1, 3) %>% 
  filter(GCF %in% intersect)

data2 <- data %>% 
  group_by(GCF) %>% 
  group_modify(~head(.x, 1))

edges <- filter(network, from %in% data$BGC & to %in% data$BGC)
nodes <- data.frame(node = unique(c(edges$from, edges$to))) %>% 
  mutate(class = ifelse(grepl("BGC", node), "MIBiG", "LiM_2024"),
         label = data2$GCF[match(node, data2$BGC)])
  
# library(tidygraph)
# library(ggraph)

graph <- tbl_graph(edges = edges, nodes = nodes)

ggraph(graph, layout = "fr") + 
  geom_edge_link(linewidth = .3, color = "grey20", linetype = "dashed") +
  geom_node_point(aes(color = class), size = 2) +
  geom_node_text(aes(label = label), size = 2) +
  scale_color_manual(values = c("#EB9256", "#58BBD0")) +
  theme_graph() +
  theme(aspect.ratio = 1)

ggsave("GCFs.vs.MIBiG.intersect.nwk.pdf", width = 7, height = 7)

# net all

files <- list.files("bigscape.out/network_files/2024-05-21_12-23-16_hybrids_glocal/", 
                    full.names = T, recursive = T, pattern = "30.network")

map(files, \(x) {
  name <- stringr::str_split_i(x, "/", 4)
  network <- read.delim(x, check.names = F) %>% 
    select(from = 1, to = 2, dist = 3)

  edges <- network
  nodes <- data.frame(node = unique(c(edges$from, edges$to))) %>% 
    mutate(class = ifelse(grepl("BGC", node), "MIBiG", "LiM_2024"),
           genome = ifelse(grepl("BGC", node), "", stringr::str_extract(node, "[A-Za-z0-9.]+"))) %>% 
    left_join(taxa %>% select(name, phylum), by = c("genome" = "name")) %>% 
    mutate(phylum = gsub("p__", "", phylum),
           phylum = gsub("_\\w$", "", phylum),
           phylum = ifelse(phylum %in% names(color_phylum), phylum, "Other"),
           label = ifelse(grepl("BGC", x = node), BGC_data$GCF[match(node, BGC_data$BGC)], ""))

  graph <- tbl_graph(edges = edges, nodes = nodes)

  color_phylum <- structure(c("#bebada","#F9C851","#8dd3c7","#FF9D9A","#fee08b","#a7c9e0","#de77ae","#fabfd2","#b8e186","#d9d9d9"),
                          names = c("Actinomycetota","Bacillota","Bacteroidota","Cyanobacteriota","Methanobacteriota", "Other","Pseudomonadota","Spirochaetota","Thermoproteota","Verrucomicrobiota"))

  ggraph(graph, layout = "circlepack") + 
    geom_edge_link(linewidth = .2, color = "grey70") +
    geom_node_point(aes(shape = class, size = class, color = phylum)) +
    scale_shape_manual(values = c(16, 17)) +
    scale_size_manual(values = c(1.2, 1.4)) +
    scale_color_manual(values = color_phylum) +
    theme_void() +
    theme(aspect.ratio = 1)

  ggsave(paste0("GCFs.nwk.", name ,".pdf"), width = 10, height = 10)
  
  })

# BGCs MIBiG GCC intersect -------------------

files <- list.files("bigscape.out/network_files/2024-05-21_12-23-16_hybrids_glocal/", 
                    full.names = T, recursive = T, pattern = "0.30_0.70.tsv")

data <- map_dfr(files, \(x) read.delim(x, check.names = F)) %>% 
  select(BGC = 1, GCC = 2, GCF = 3) %>% 
  mutate(GCC = paste0("GCC", sprintf("%05d", GCC)),
         GCF = paste0("GCF", sprintf("%05d", GCF)),
         class = ifelse(grepl("BGC", BGC), "MIBiG", "LiM_2024"))

data$GCC %>% unique() %>% length

plot_data <- data %>% 
  group_by(GCC) %>% 
  group_modify(~.x %>%
                 group_by(class) %>% 
                 summarise(n = n())) 

plot_data %>% 
  spread(key = "class", value = "n", fill = 0) %>% 
  column_to_rownames(var = "GCC") %>%
  apply(2, \(x) ifelse(x>0, rownames(.), NA)) %>% 
  data.frame %>% 
  map(\(x) as.character(na.omit(x))) %>% 
  ggvenn::ggvenn(., fill_color = c("#f1b389", "#8ad0de"), stroke_linetype = "dashed", 
                 auto_scale = T, fill_alpha = .9)
ggsave("GCCs.vs.MIBiG.venn.pdf", width = 5, height = 4)

# nwk
plot_data <- data %>% 
  group_by(GCC) %>% 
  group_modify(~.x %>%
                 group_by(class) %>% 
                 summarise(n = n())) 

intersect <- plot_data %>% 
  spread(key = "class", value = "n", fill = 0) %>% 
  column_to_rownames(var = "GCC") %>%
  filter(MIBiG > 0 & LiM_2024 > 0) %>% 
  rownames()

files <- list.files("bigscape.out/network_files/2024-05-21_12-23-16_hybrids_glocal/", 
                    full.names = T, recursive = T, pattern = "30.network")

network <- map_dfr(files, \(x) read.delim(x, check.names = F)) %>% 
  select(from = 1, to = 2, dist = 3)

data <- data %>% 
  select(1, 2) %>% 
  filter(GCC %in% intersect)

data2 <- data %>% 
  group_by(GCC) %>% 
  group_modify(~head(.x, 1))

edges <- filter(network, from %in% data$BGC & to %in% data$BGC)
nodes <- data.frame(node = unique(c(edges$from, edges$to))) %>% 
  mutate(class = ifelse(grepl("BGC", node), "MIBiG", "LiM_2024"),
         GCC = data2$GCC[match(node, data2$BGC)])
  
graph <- tbl_graph(edges = edges, nodes = nodes)

ggraph(graph, layout = "graphopt") + 
  geom_edge_link(linewidth = .1, color = "grey20") +
  geom_node_point(aes(color = class), size = .3, position = position_jitter()) +
  scale_color_manual(values = c("#EB9256", "#58BBD0")) +
  theme_graph() +
  theme(aspect.ratio = 1)
ggsave("GCCs.vs.MIBiG.intersect.nwk.pdf", width = 7, height = 7)

# BGCs filter the longest BGC per GCF per species ------------------------------
# refer to "Biosynthetic potential of the global ocean microbiome" published on Nature.
# To prevent sampling biases in quantitative analysis (taxonomic and functional 
# compositions of GCCs/GCFs, GCF and GCC distances to reference databases as well
# as GCF metagenomic abundances), the 73,864 BGCs were further dereplicated by 
# retaining only the longest BGC per GCF per species, 
# resulting in a total of 30,182 BGCs.
# 有些BGC可能归为2中GCF甚至是GCC ，unique的BGCs是15919，包括归为多种GCF的BGC行数是16323

# file_path <- "bigscape.out/network_files/2024-05-21_12-23-16_hybrids_glocal/"
# files <- paste0(file_path, dir(file_path, recursive = T, pattern = "0.30_0.70.tsv"))
# data <- map_dfr(files, \(x) read.delim(x, check.names = F)) %>%
#   select(BGC = 1, GCC = 2, GCF = 3) %>%
#   filter(!grepl("BGC", BGC)) %>%
#   mutate(GCC = paste0("GCC", sprintf("%05d", GCC)),
#          GCF = paste0("GCF", sprintf("%05d", GCF)),)
# 
# data$GCC %>% unique() %>% length()
# data$GCF %>% unique() %>% length()
# 
# data2 <- read.delim("rawdata/Network_Annotations_Full.tsv") %>%
#   select(BGC, accession = Accession.ID, description = Description, product_prediction = "Product.Prediction", BiGSCAPE_class = "BiG.SCAPE.class")
# 
# data3 <- read.delim("../prok_genome/data.strain.genome.metadata.txt") %>%
#   select(name, domain, phylum, family, genus, species, genomospecies)
# 
# data4 <- read.delim("rawdata/BGCs.gbk.len", header = F, col.names = c("name", "gbk", "length")) %>%
#   mutate(BGC = gsub(".gbk", "", gbk)) %>%
#   select(BGC, length)
# 
# data %>%
#   mutate(genome = gsub("_.*", "", BGC)) %>%
#   left_join(data2, by = "BGC") %>%
#   left_join(data4, by = "BGC") %>%
#   left_join(data3, by = c("genome" = "name")) %>%
#   relocate(genomospecies, .after = "species")

metadata <- read.delim("BGCs.metadata.tsv")
data <- metadata %>% 
  group_by(genomospecies, GCF) %>% 
  group_modify(~slice_max(.x, order_by = length, with_ties = F))
write.table(data, "BGCs.metadata.select_BGCs_per_species_per_GCFs.tsv", sep = "\t", row.names = F, quote = F)

data2 <- metadata %>% 
  group_by(genomospecies, GCF) %>% 
  summarise(n = n())
write.table(data2, "BGCs.metadata.select_BGCs_per_species_per_GCFs.count.tsv", sep = "\t", row.names = F, quote = F)

data %>% 
  ungroup %>% 
  select(BGC) %>% 
  unique() %>% 
  write.table("BGCs.metadata.select_BGCs_per_species_per_GCFs.list", col.names = F, row.names = F, quote = F)

# BGCs novo -----------------------------------------------------------------------------
# 有些BGC可能归为2中GCF甚至是GCC ，unique的BGCs是15919，包括归为多种GCF的BGC行数是16323

BGCs_rep <- read.delim("BGCs.metadata.select_BGCs_per_species_per_GCFs.tsv")

# BGC
BGCs_unmap <- read.delim("rawdata/bigslice.out.metadata.d_m900", header = F) %>% 
  select(x = 5) %>% 
  mutate(BGC = stringr::str_split_i(x, "/", 2))

# 3436/7371  # novo / total
# nrow(BGCs_unmap) / length(unique(BGCs_rep$BGC))

# 不同类型BGC新型比例
plot_data <- BGCs_rep %>% 
  select(BGC, BiGSCAPE_class) %>% 
  unique.data.frame %>% 
  mutate(novo = ifelse(BGC %in% BGCs_unmap$BGC, "novo", "known")) %>% 
  group_by(BiGSCAPE_class, novo) %>% 
  summarise(n = n()) %>% 
  spread(key = "novo", value = "n", fill = 0) %>% 
  mutate(novo_prec = novo / (known + novo)) %>% 
  data.frame

plot_data <- plot_data %>% 
  mutate(sum = known + novo) %>% 
  arrange(sum) %>% 
  mutate(BiGSCAPE_class = factor(BiGSCAPE_class, BiGSCAPE_class))

colors <- c("#82ccdc","#f0ad80","#b4b3d8","#77b0c5","#c2b75f","#87CCBA","#FEDFB2","#f08b85")
names(colors) <- c("RiPPs","Others","NRPS","Terpene","PKSother","PKSI","PKS-NRP_Hybrids","Saccharides")

p1 <- ggbarplot(plot_data, x = "BiGSCAPE_class", y = "sum", fill = "BiGSCAPE_class", 
                xlab = "", sorting = "none", rotate = T, ylab = "Number of BGCs", 
                label = T, lab.vjust = .5, lab.hjust = -.1, palette = colors, width = .61) +
  scale_y_continuous(position = "right") +
  theme(aspect.ratio = 4/3)

plot_label <- plot_data %>% 
  mutate(BiGSCAPE_class = factor(BiGSCAPE_class, plot_data$BiGSCAPE_class)) %>% 
  mutate(label = paste0("Novel BGCs:", round(novo_prec*100, 2), "%"))

p2 <- ggplot(plot_data, aes(x = BiGSCAPE_class, y = 100, fill = BiGSCAPE_class)) +
  geom_bar(stat = "identity", color = "#000000", width = .61, alpha = .5) +
  geom_bar(aes(x = BiGSCAPE_class, y = novo_prec * 100, fill = BiGSCAPE_class), plot_data,
           stat = "identity", inherit.aes = F, width = .61, color = "#000000") +
  geom_label(aes(x = BiGSCAPE_class, y = novo_prec * 100 - 1, label = label, hjust = 1), 
             plot_label, size = 3, color = "white", fill = "transparent") +
  scale_fill_manual(values = colors) +
  scale_y_continuous(position = "right") +
  labs(x = "", y = "Novo") +
  coord_flip() +
  theme_pubr() +
  theme(aspect.ratio = 4/3)

cowplot::plot_grid(p2, p1, align = "v")

ggsave("BGCs.novo.prec.pdf", width = 10, height = 6)

# 不同物种来源BGC新型比例

BGCs_rep %>% 
  select(BGC, phylum) %>% 
  mutate(phylum = gsub("_\\w$", "", phylum)) %>% 
  unique.data.frame %>% 
  mutate(novo = ifelse(BGC %in% BGCs_unmap$BGC, "novo", "known")) %>% 
  group_by(phylum, novo) %>% 
  summarise(n = n()) %>% 
  spread(key = "novo", value = "n", fill = 0) %>% 
  mutate(novo_prec = novo / (known + novo)) %>% 
  data.frame %>% 
  write.table("BGCs.novo.prec.phylum.txt", sep = "\t", quote = F, row.names = F)

# GCFs novo --------------

data <- BGCs_rep %>% 
  select(GCF, BGC) %>%
  mutate(novo = ifelse(BGC %in% BGCs_unmap$BGC, "unmapped", "mapped")) %>% 
  group_by(GCF) %>% 
  group_modify(~table(.x$novo) %>% data.frame)

data2 <- data %>% 
  spread(key = "Var1", value = "Freq", fill = 0) %>% 
  mutate(prec = mapped / (mapped + unmapped)) %>% 
  filter(prec < .2)

# 2275 / 5299 # novo / total
# nrow(dat2) / length(unique(BGCs_rep$GCF))

# 不同类型GCF新型比例
plot_data <- BGCs_rep %>% 
  select(GCF, BiGSCAPE_class) %>% 
  distinct() %>% 
  mutate(novo = ifelse(GCF %in% data2$GCF, "novo", "known")) %>% 
  group_by(BiGSCAPE_class, novo) %>% 
  summarise(n = n()) %>% 
  spread(key = "novo", value = "n", fill = 0) %>% 
  mutate(novo_prec = novo / (known + novo)) %>% 
  data.frame

plot_data <- plot_data %>% 
  mutate(sum = known + novo) %>% 
  arrange(sum) %>% 
  mutate(BiGSCAPE_class = factor(BiGSCAPE_class, BiGSCAPE_class))

p1 <- ggbarplot(plot_data, x = "BiGSCAPE_class", y = "sum", fill = "BiGSCAPE_class", 
                xlab = "", sorting = "none", rotate = T, ylab = "Number of GCFs", 
                label = T, lab.vjust = .5, lab.hjust = -.1, palette = colors, width = .61) +
  scale_y_continuous(position = "right") +
  theme(aspect.ratio = 4/3)

plot_label <- plot_data %>% 
  mutate(BiGSCAPE_class = factor(BiGSCAPE_class, plot_data$BiGSCAPE_class)) %>% 
  mutate(label = paste0("Novel GCFs:", round(novo_prec*100, 2), "%"))

p2 <- ggplot(plot_data, aes(x = BiGSCAPE_class, y = 100, fill = BiGSCAPE_class)) +
  geom_bar(stat = "identity", color = "#000000", width = .61, alpha = .5) +
  geom_bar(aes(x = BiGSCAPE_class, y = novo_prec * 100, fill = BiGSCAPE_class), plot_data,
           stat = "identity", inherit.aes = F, width = .61, color = "#000000") +
  geom_label(aes(x = BiGSCAPE_class, y = novo_prec * 100 - 1, label = label, hjust = 1), 
             plot_label, size = 3, color = "white", fill = "transparent") +
  scale_fill_manual(values = colors) +
  scale_y_continuous(position = "right") +
  labs(x = "", y = "Novo") +
  coord_flip() +
  theme_pubr() +
  theme(aspect.ratio = 4/3)

cowplot::plot_grid(p2, p1, align = "v")

ggsave("GCFs.novo.prec.pdf", width = 10, height = 6)

# 不同物种来源GCF新型比例
BGCs_rep %>% 
  select(GCF, phylum) %>% 
  mutate(phylum = gsub("_\\w$", "", phylum)) %>% 
  distinct() %>% 
  mutate(novo = ifelse(GCF %in% data2$GCF, "novo", "known")) %>% 
  group_by(phylum, novo) %>% 
  summarise(n = n()) %>% 
  spread(key = "novo", value = "n", fill = 0) %>% 
  mutate(novo_prec = novo / (known + novo)) %>% 
  data.frame %>% 
  write.table("GCFs.novo.prec.phylum.txt", sep = "\t", quote = F, row.names = F)

# GCCs novo ---------------------

data <- BGCs_rep %>% 
  select(GCC, BGC) %>%
  mutate(novo = ifelse(BGC %in% BGCs_unmap$BGC, "unmapped", "mapped")) %>% 
  group_by(GCC) %>% 
  group_modify(~table(.x$novo) %>% data.frame)

data2 <- data %>% 
  spread(key = "Var1", value = "Freq", fill = 0) %>% 
  mutate(prec = mapped / (mapped + unmapped)) %>% 
  filter(prec < .4)

# 265 / 847 # novo / total
# nrow(dat2) / length(unique(BGCs_rep$GCC))

# 不同类型GCC新型比例
plot_data <- BGCs_rep %>% 
  select(GCC, BiGSCAPE_class) %>% 
  distinct() %>% 
  mutate(novo = ifelse(GCC %in% data2$GCC, "novo", "known")) %>% 
  group_by(BiGSCAPE_class, novo) %>% 
  summarise(n = n()) %>% 
  spread(key = "novo", value = "n", fill = 0) %>% 
  mutate(novo_prec = novo / (known + novo)) %>% 
  data.frame

plot_data <- plot_data %>% 
  mutate(sum = known + novo) %>% 
  arrange(sum) %>% 
  mutate(BiGSCAPE_class = factor(BiGSCAPE_class, BiGSCAPE_class))

p1 <- ggbarplot(plot_data, x = "BiGSCAPE_class", y = "sum", fill = "BiGSCAPE_class", 
                xlab = "", sorting = "none", rotate = T, ylab = "Number of GCCs", 
                label = T, lab.vjust = .5, lab.hjust = -.1, palette = colors, width = .61) +
  scale_y_continuous(position = "right") +
  theme(aspect.ratio = 4/3)

plot_label <- plot_data %>% 
  mutate(BiGSCAPE_class = factor(BiGSCAPE_class, plot_data$BiGSCAPE_class)) %>% 
  mutate(label = paste0("Novel GCCs:", round(novo_prec*100, 2), "%"))

p2 <- ggplot(plot_data, aes(x = BiGSCAPE_class, y = 100, fill = BiGSCAPE_class)) +
  geom_bar(stat = "identity", color = "#000000", width = .61, alpha = .5) +
  geom_bar(aes(x = BiGSCAPE_class, y = novo_prec * 100, fill = BiGSCAPE_class), plot_data,
           stat = "identity", inherit.aes = F, width = .61, color = "#000000") +
  geom_label(aes(x = BiGSCAPE_class, y = novo_prec * 100 - 1, label = label, hjust = 1), 
             plot_label, size = 3, color = "white", fill = "transparent") +
  scale_fill_manual(values = colors) +
  scale_y_continuous(position = "right") +
  labs(x = "", y = "Novo") +
  coord_flip() +
  theme_pubr() +
  theme(aspect.ratio = 4/3)

cowplot::plot_grid(p2, p1, align = "v")

ggsave("GCCs.novo.prec.pdf", width = 10, height = 6)

# 不同物种来源GCC新型比例
BGCs_rep %>% 
  select(GCC, phylum) %>% 
  mutate(phylum = gsub("_\\w$", "", phylum)) %>% 
  distinct() %>% 
  mutate(novo = ifelse(GCC %in% data2$GCC, "novo", "known")) %>% 
  group_by(phylum, novo) %>% 
  summarise(n = n()) %>% 
  spread(key = "novo", value = "n", fill = 0) %>% 
  mutate(novo_prec = novo / (known + novo)) %>% 
  data.frame %>% 
  write.table("GCCs.novo.prec.phylum.txt", sep = "\t", quote = F, row.names = F)

# profile BGCs -----------------------------------------------------------------------------

profile <- read.delim("rawdata/profile.rc3") %>% 
  mutate(name = stringr::str_extract(name, "[\\w.]+")) %>% 
  aggregate(. ~ name, ., sum) %>% 
  column_to_rownames(var = "name")

lib_size <- read.delim("rawdata/profile.lib_size", col.names = c("name", "lib_size"), header = F) %>% 
  filter(name %in% colnames(profile)) %>% 
  mutate(name = factor(name, colnames(profile))) %>% 
  arrange(name)

rela_ab <- apply(profile, 1, \(x) x / lib_size$lib_size * 100) %>% 
  t() %>% 
  data.frame

saveRDS(rela_ab, "rawdata/profile.rela.ab.rds")

# div for breed --------------------

source("/code/R_func/diversity.R")
group <- read.delim("../profile/sample_group")
profile <- readRDS("rawdata/profile.rela.ab.rds")
metdata <- read.delim("BGCs.metadata.tsv")

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

data <- calcu_alpha(profile_x, "richness") %>% 
  left_join(group_x, by = "sample")

var <- data %>% 
  filter(grepl("PN5R", sample)) %>% 
  pull(2) %>% 
  median()

p1 <- ggboxplot(data, x = "group", y = "val", xlab = "", ylab = "Richness Index", 
                fill = "class", legend = "none", rotate = T, color = "black",
                outlier.shape = NA, palette = c("#56B4E9", "#E69F00")) +
  stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", hide.ns = T, color = "red") +
  geom_hline(yintercept = var, lty = "dashed", linewidth = .4) +
  theme(aspect.ratio = 2)

data <- calcu_alpha(profile_x, "shannon") %>% 
  left_join(group_x, by = "sample")

var <- data %>% 
  filter(grepl("PN5R", sample)) %>% 
  pull(2) %>% 
  median()

p2 <- ggboxplot(data, x = "group", y = "val", xlab = "", ylab = "Shannon Index", 
                fill = "class", legend = "none", rotate = T, color = "black",
                outlier.shape = NA, palette = c("#56B4E9", "#E69F00")) +
  stat_compare_means(ref.group = "JL-DLLW", label = "p.signif", hide.ns = T, color = "red") +
  geom_hline(yintercept = var, lty = "dashed", linewidth = .4) +
  theme(aspect.ratio = 2)

cowplot::plot_grid(p1, p2, align = "v")
ggsave("alpha_div.for_breed.pdf", height = 5, width = 8)

# div for breed meta --------------------------------------------------------------------------------------------------------

source("/code/R_func/calcu_metafor.R")
source("/code/R_func/diversity.R")
source("/code/R_func/difference_analysis.R")

metdata <- read.delim("BGCs.metadata.tsv")

data <- calcu_alpha(profile_x, "richness") %>% 
  left_join(group_x, by = "sample")

vec_ref <- data[data$group == "JL-DLLW", "val"]

meta_out <- data %>% 
  filter(group != "JL-DLLW") %>% 
  select(val, group) %>% 
  group_by(group) %>% 
  group_modify(~data.frame(mean = mean(.x$val), sd = sd(.x$val), n = length(.x$val))) %>%
  ungroup %>% 
  rename(d_mean = mean, d_sd = sd, d_n = n) %>% 
  add_column(c_mean = mean(vec_ref), c_sd = sd(vec_ref), c_n = length(vec_ref)) %>% 
  rename(name = group) %>% 
  metafor_fit.1() %>% 
  add_column(index = "richness")

data <- calcu_alpha(profile_x, "shannon") %>% 
  left_join(group_x, by = "sample")

vec_ref <- data[data$group == "JL-DLLW", "val"]

meta_out <- rbind(
  meta_out, 
  data %>% 
    filter(group != "JL-DLLW") %>% 
    select(val, group) %>% 
    group_by(group) %>% 
    group_modify(~data.frame(mean = mean(.x$val), sd = sd(.x$val), n = length(.x$val))) %>%
    ungroup %>% 
    rename(d_mean = mean, d_sd = sd, d_n = n) %>% 
    add_column(c_mean = mean(vec_ref), c_sd = sd(vec_ref), c_n = length(vec_ref)) %>% 
    rename(name = group) %>% 
    metafor_fit.1() %>% 
    add_column(index = "shannon")
)

write.table(meta_out, "alpha_div.for_breed.meta.txt", sep = "\t", quote = F, row.names = F)

meta_out <- meta_out %>% 
  mutate(padj = p.adjust(pval, "BH")) %>% 
  add_plab(by = "padj")

ggplot() +
  geom_vline(xintercept = 0, color = "#000000", lwd = .4, lty = 2) +
  geom_point(aes(x = estimate, y = index), 
             meta_out %>% select(index, estimate) %>% unique(), 
             size = 2.5, shape = 15, inherit.aes = F) + 
  geom_errorbar(aes(y = index, xmin = ci_lb, xmax = ci_ub),
                meta_out %>% select(index, ci_lb, ci_ub) %>% unique(),
                width = .2, lwd = .4, inherit.aes = F) +
  geom_point(aes(x = yi, y = index, color = name),
             meta_out, 
             size = 2, inherit.aes = F, show.legend = F) +
  geom_text(aes(x = min(yi), y = index, label =  plab), 
            meta_out, color = "red") +
  labs(y = "", x = "Standradized Mean Difference (Random Effect Model)") +
  scale_color_viridis_d(begin = .4) +
  theme_classic2() +
  theme(axis.text = element_text(color = "black", size = 8),
        axis.title = element_text(size = 8),
        axis.line = element_line(linewidth = .4, color = "#000000"),
        axis.ticks = element_line(linewidth = .4, color = "#000000"))

ggsave("alpha_div.for_breed.meta.pdf", width = 8, height = 2)

# div for region ----------------------

source("F:/code/R_func/diversity.R")
source("/code/R_func/calcu_metafor.R")
group <- read.delim("../profile/sample_group")
profile <- readRDS("rawdata/profile.rela.ab.rds")

region_order <- c("Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Rectum")

group_x <- group %>% 
  filter(breed2 %in% c("JL-SLBP", "JL-LXP", "JL-BSWP", "JL-DLLW")) %>%
  select(sample, breed = breed2, region)

profile_x <- profile %>% 
  select(any_of(group_x$sample)) %>% 
  filter(rowSums(.) != 0)

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

cowplot::plot_grid(p1, p2, align = "v")
ggsave("alpha_div.for_region.pdf", height = 6, width = 8)

# PCoA for breed ----------------------------

library(patchwork)
source("/code/R_func/plot_PCoA.R")
source("/code/R_func/difference_analysis.R")
group <- read.delim("../profile/sample_group")
profile <- readRDS("rawdata/profile.rela.ab.rds")
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
  write.table("beta_div.bray.for_breed.disper.txt", sep = "\t", quote = F, row.names = F)

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

ggsave("beta_div.bray.for_breed.manual.pdf", height = 12, width = 14)

# PCoA for region facet by breed -----------------

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
ggsave("beta_div.bray.for_region.facet_by_breed.pdf", height = 3, width = 12)

cowplot::plot_grid(plotlist = plot_list2, nrow = 1)
ggsave("beta_div.bray.for_region.facet_by_breed.boxplot.pdf", height = 4, width = 15)

write.table(adonis2_data, "beta_div.bray.for_region.facet_by_breed.adonis.txt", 
            sep = "\t", quote = F, row.names = F)

adonis2_data %>% 
  mutate(name = factor(name, breed_order),
         plab = ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))) %>% 
  ggbarplot(x = "name", y = "r2adj", label = .$plab, fill = "name", palette = colors_breed,
            legend = "none", color = "black", x.text.angle = 90, lab.col = "red")
ggsave("beta_div.bray.for_region.facet_by_breed.adonis.barplot.pdf", height = 5, width = 5)

disper_data %>% 
  write.table("beta_div.bray.for_region.facet_by_breed.disper.txt", sep = "\t", row.names = F, quote = F)
disper_data2 %>% 
  write.table("beta_div.bray.for_region.facet_by_breed.disper2.txt", sep = "\t", row.names = F, quote = F)

# GCF GCC profile ----------------

metadata <- read.delim("BGCs.metadata.tsv")

profile <- readRDS("rawdata/profile.rela.ab.rds") %>%
  mutate(GCF = metadata$GCF[match(rownames(.), metadata$BGC)]) %>% 
  aggregate(. ~ GCF, ., sum) %>% 
  column_to_rownames(var = "GCF")

saveRDS(profile, "rawdata/profile.GCF.rela.ab.rds")

profile <- readRDS("rawdata/profile.rela.ab.rds") %>%
  mutate(GCC = metadata$GCC[match(rownames(.), metadata$BGC)]) %>% 
  aggregate(. ~ GCC, ., sum) %>% 
  column_to_rownames(var = "GCC")

saveRDS(profile, "rawdata/profile.GCC.rela.ab.rds")

profile <- readRDS("rawdata/profile.rela.ab.rds") %>%
  mutate(BiGSCAPE_class = metadata$BiGSCAPE_class[match(rownames(.), metadata$BGC)]) %>% 
  aggregate(. ~ BiGSCAPE_class, ., sum) %>% 
  column_to_rownames(var = "BiGSCAPE_class")

saveRDS(profile, "rawdata/profile.BiGSCAPE_class.rela.ab.rds")


  
