library(here)
library(ape)
library(ggtree)
library(dplyr)

print("Reading tree.")
tree <- read.tree(here("../radseq_pruned_csubst.tre"))

p <- ggtree(tree) +
  geom_tiplab(hjust = -0.025, size = 2, fontface = "italic", align = TRUE, offset = 0.001) +
  theme(legend.position = "none") +
  hexpand(0.15, direction = 1) +
  xlim_tree(c(-0.0025, 0.035))

edge <- data.frame(tree$edge, edge_num=1:nrow(tree$edge))
colnames(edge) <- c("parent", "node", "edge_num")
edge <- edge %>% filter(edge_num %in% c(160, 174, 228))

p %<+% edge + geom_label(aes(x = branch, label = edge_num), size = 2)
