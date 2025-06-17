library(dplyr)
here::i_am("5. unroot_tree.R")
tree <- ape::read.tree("pruned.tre") # read tree
tree$tip.label <- tree$tip.label %>%
  stringr::str_remove_all("'") %>%
  stringr::str_remove_all("_E634") %>%
  stringr::str_remove_all("_E492") %>%
  stringr::str_remove_all("_E489") %>%
  stringr::str_remove_all("_Diacalpe|_Acrophorus|_Nothoperanema|_Peranema|_Dryopsis|'")

sessa_accessions <- readr::read_csv("sessa_tree_accessions.csv") %>% # in case I missed labels
  rename(accession = rbcL) %>%
  select(accession, ID) %>%
  filter(accession != "—", ID != "—")


rbcL_names <- readr::read_csv("metadata.csv") %>%
  left_join(sessa_accessions) %>%
  mutate(species = case_when(!stringr::str_detect(species, ID) & !is.na(ID) ~ stringr::str_c(species, ID, sep = " "),
                          .default = species)) %>%
  mutate(TreeName = stringr::str_replace_all(species, " ", "_"),
         TreeName = stringr::str_replace_all(TreeName, "-", "_"),
         NewName = case_when(
           stringr::str_detect(species, "\\sE[:digit:]{3}$") ~ stringr::str_c(stringr::str_sub(stringr::word(species, 2), 1, 7), 
                                                                       stringr::str_sub(stringr::word(species, 3), 2, 4)),
           .default = stringr::str_sub(stringr::word(species, 2), 1, 10))) %>%
  mutate(TreeName = stringr::str_replace_all(TreeName, "azorica", "intermedia_subsp._azorica"),
         TreeName = stringr::str_replace_all(TreeName, "coreano_montana", "coreanomontana"),
         TreeName = stringr::str_replace_all(TreeName, "fuscoatra", "fusco_atra"),
         TreeName = stringr::str_remove(TreeName, "_E505"),
         TreeName = stringr::str_replace_all(TreeName, "monticola", "goldiana_ssp._monticola"),
         TreeName = stringr::str_replace_all(TreeName, "maderensis", "intermedia_subsp._maderensis"),
         TreeName = stringr::str_replace_all(TreeName, "sabaei", "sabae")) %>%
  filter(!stringr::str_detect(TreeName, "var._soripes")) %>%
  select(TreeName, NewName)

rbcL_taxa <- rbcL_names %>% pull(TreeName) %>% sort()

pruned_tree <- ape::keep.tip(tree, rbcL_taxa) %>% 
  treeio::rename_taxa(tree = ., data = rbcL_names, key = TreeName, value = NewName)

ape::write.tree(pruned_tree, file = "dryopteris_csubst.tre")
pruned_tree$edge.length <- NULL # remove branch lengths for 
ape::write.tree(ape::unroot(pruned_tree), file = "dryopteris_paml.tre")

# 
# temp <- ape::read.tree("dryopteris_csubst.tre")$tip.label
# temp[duplicated(temp)]
