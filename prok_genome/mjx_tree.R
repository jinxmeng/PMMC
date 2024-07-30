# Jinxin Meng, 20240125, 20240707 --------------------------

pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr, forcats)
setwd("F:/proj/proj_2024/20240425_metagenome_Black_pigs_Limh/prok_genome/")

# itol ---------------------------

# 基因组流行率

source("/code/R_func/profile_process.R")

group <- read.delim("../profile/sample_group")
profile <- readRDS("../profile/genomospecies.tpm.b50.rds")

group_x <- group %>%
  filter(group2 != "")

profile_x <- profile %>% 
  select(any_of(group_x$sample))

data <- profile_prevalence(profile_x, by_group = F)
write.table(data, "tree/prok.occurrence.tsv", sep = "\t", row.names = F, quote = F)

data <- profile_x %>% 
  rowMeans() %>% 
  data.frame(value = .) %>% 
  rownames_to_column("name")
write.table(data, "tree/prok.abundance.tsv", sep = "\t", row.names = F, quote = F)

# ggtree 物种进化树 -------------------------------
# 
# library(ggtree)
# library(ggtreeExtra)
# library(tidytree)
# library(ggnewscale)
# library(treeio)
# 
# source("/code/R_func/profile_process.R")
# source("/code/R_func/utilities.R")
# 
# metadata <- read.delim("data.species.genome.metadata.txt")
# 
# tr <- read.tree("genomospecies.tre") %>%
#   root(outgroup = "PN6RQ.bin.90", edgelabel = T)
# 
# tr_data <- as_tibble(tr)
# 
# tr_bar_phylum <- select(metadata, genome, phylum) %>%
#   mutate(phylum = fct_lump_n(phylum, 10, other_level = "p__Other"))
# tr_phylum_colors <- tr_bar_phylum$phylum %>%
#   table %>%
#   sort(decreasing = T) %>%
#   names() %>%
#   discard(~ . == "p__Other") %>%
#   append("p__Other") %>%
#   structure(c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#ACDEF3",
#               "#F9C6B3","#F5EAF2","#F9C851","#D3EDE4", "#d9d9d9"), names = .)
# 
# tr_bar_source <- select(metadata, genome, source)
# tr_source_colors <- structure(c("#ee7e77","#ffffff"), names = c("isolate", "metagenome-assembly"))
# 
# tr_bar_expand <- select(metadata, genome, expand)
# tr_expand_colors <- structure(c("#66c2a5","#ffffff"), names = c("Yes", "No"))
# 
# tr_bar_novo <- select(metadata, genome, novo)
# tr_novo_colors <- structure(c("#e69f00","#ffffff"), names = c("Yes", "No"))
# 
# tr_bar_count <- select(metadata, genome, n = strains)
# 
# p <- ggtree(tr, layout = "fan", lwd = .4, color = "grey60") +
#   geom_rootedge(rootedge = .06) +
#   geom_fruit(aes(y = genome, fill = phylum),
#              tr_bar_phylum, 
#              geom_tile, 
#              pwidth = .09, offset = .15, inherit.aes = F) +
#   scale_fill_manual(values = tr_phylum_colors) +
#   new_scale_fill() +
#   geom_fruit(aes(y = genome, fill = source),
#              tr_bar_source,
#              geom_tile, 
#              pwidth = .09, offset = .15, inherit.aes = F) +
#   scale_fill_manual(values = tr_source_colors) +
#   new_scale_fill() +
#   geom_fruit(aes(y = genome, fill = expand),
#              tr_bar_expand,
#              geom_tile,
#              pwidth = .09, offset = .15, inherit.aes = F) +
#   scale_fill_manual(values = tr_expand_colors) +
#   new_scale_fill() +
#   geom_fruit(aes(y = genome, fill = novo),
#              tr_bar_novo,
#              geom_tile,
#              pwidth = .09, offset = .15, inherit.aes = F) +
#   scale_fill_manual(values = tr_novo_colors) +
#   new_scale_fill() +
#   geom_fruit(aes(y = genome, x = n),
#              tr_bar_count, 
#              geom_bar, 
#              stat = "identity", orientation = "y", pwidth = .4, offset = .15, 
#              inherit.aes = F, grid.params = list(color = "#56b4e9", size = .6))
# open_tree(p, angle = 15)
# ggsave("spcies.ggtree.v2.pdf", width = 20, height = 20)  
