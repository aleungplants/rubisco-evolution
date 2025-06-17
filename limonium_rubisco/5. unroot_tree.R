library(dplyr)
here::i_am("5. unroot_tree.R")
tree <- ape::read.tree("koutroumpa_tree.tre") # read tree
tree$tip.label <- tree$tip.label %>%
  stringr::str_replace_all("recurvum_subsp_humile", "rhumile") %>%
  stringr::str_replace_all("vulgare_eduardi_diasii", "diasii") %>%
  stringr::str_replace_all("_", " ") %>%
  stringr::word(-1) %>%
  stringr::str_sub(1, 10) # change tip labels to match short phylip names
rbcL_taxa <- readr::read_csv("accessions.csv") %>%
  filter(!is.na(accession)) %>%
  mutate(species = stringr::str_replace(species, "vulgare \\(eduardi-diasii", "diasii")) %>%
  mutate(species = stringr::str_replace(species, "recurvum C.E.Salmon subsp. humile", "rhumile")) %>%
  mutate(species = stringr::word(species, 2, sep = "(?=[A-Z()])")) %>%
  pull(species) %>%
  stringr::str_replace_all(" cf ", " ") %>%
  stringr::str_replace_all("-", " ") %>%
  stringr::str_remove_all(" $") %>%
  stringr::str_replace_all("barba", "jovibarba") %>%
  stringr::word(-1) %>%
  stringr::str_sub(1, 10) %>%
  stringr::str_remove_all("\\.")
taxa_keep <- intersect(tree$tip.label, rbcL_taxa)
tree <- ape::keep.tip(tree, taxa_keep)
ape::write.tree(tree, file = "koutroumpa_tree_csubst.tre")
tree$edge.length <- NULL # remove branch lengths for 
ape::write.tree(ape::unroot(tree), file = "koutroumpa_tree_paml.tre")
