library(dplyr)
here::i_am("5. unroot_tree.R")
tree <- ape::read.tree("Pinus_time.tre") # read tree
tree$tip.label <- tree$tip.label %>%
  stringr::str_replace_all("_", " ") %>%
  stringr::word(2) %>%
  stringr::str_sub(1, 10) # change tip labels to match short phylip names
rbcL_taxa <- readr::read_csv("accessions.csv") %>%
  pull(Species) %>%
  stringr::word(2) %>%
  stringr::str_sub(1, 10)
tree <- ape::keep.tip(tree, rbcL_taxa)
ape::write.tree(tree, file = "Pinus_time_csubst.tre")
tree$edge.length <- NULL # remove branch lengths for 
ape::write.tree(ape::unroot(tree), file = "Pinus_time_paml.tre")
