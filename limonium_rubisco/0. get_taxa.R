library(dplyr)

here::i_am("0. get_taxa.R")
limonium_tree <- ape::read.tree(here::here("koutroumpa_tree.tre")) # Koutroumpa et al 2021 front plant sci
limonium_species <- limonium_tree$tip.label %>%
  stringr::str_replace_all("_", " ") %>%
  as_tibble()

readr::write_csv(limonium_species, here::here("koutroumpa_species.csv"), col_names = FALSE)
