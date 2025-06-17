library(dplyr)

here::i_am("0. get_taxa.R")
pinus_tree <- ape::read.tree(here::here("Pinus_time.tre")) # Jin et al. 2021 PNAS Fig. 2 tree
pinus_species <- pinus_tree$tip.label %>%
  stringr::str_replace_all("_", " ") %>%
  stringr::word(1, 2) %>%
  as_tibble()

readr::write_csv(pinus_species, here::here("jin_species.csv"), col_names = FALSE)
