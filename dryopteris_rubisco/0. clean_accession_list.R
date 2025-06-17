library(dplyr)
here::i_am("0. clean_accession_list.R")
tree_taxa <- ape::read.tree(here::here("pruned.tre"))$tip.label %>%
  stringr::str_replace_all("filix_mas", "filix-mas") %>%
  stringr::str_replace_all("fusco_atra", "fuscoatra") %>%
  stringr::str_replace_all("_", " ") %>%
  stringr::str_remove_all(" Diacalpe| Acrophorus| Nothoperanema| Peranema| Dryopsis|'")
taxa_labels <- tree_taxa %>%
  as_tibble() %>%
  rename(Taxon = value) %>%
  filter(!stringr::str_detect(Taxon, "Out")) %>%
  mutate(Subsp_var = stringr::str_extract(Taxon, "(var.|subsp.)\\s[:alnum:]+"),
         Taxon = stringr::word(Taxon, 1, 2))
readr::write_csv(taxa_labels, here::here("sessa_tree_labels.csv"))
