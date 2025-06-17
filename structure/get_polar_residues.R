library(dplyr)

here::i_am("structure/get_polar_residues.R")

h_bond_v <- readxl::read_xlsx(here::here("analyze/h_bonds.xlsx")) %>%
  mutate(Genus = "Viburnum")
h_bond_d <- readxl::read_xlsx(here::here("dryopteris_rubisco/h_bonds.xlsx")) %>%
  mutate(Genus = "Dryopteris")
h_bond_l <- readxl::read_xlsx(here::here("limonium_rubisco/h_bonds.xlsx")) %>%
  mutate(Genus = "Limonium")
h_bond_p <- readxl::read_xlsx(here::here("pine_rubisco/h_bonds.xlsx")) %>%
  mutate(Genus = "Pinus")

h_bond <- bind_rows(h_bond_v, h_bond_d, h_bond_p) %>%
  arrange(Site) %>%
  mutate(Substitution = paste0(AA_Ancestral, Site, "→", AA),
         BondableChange = case_when(AA_HBondable == FALSE & AA_A_HBondable == TRUE ~ "Removed", 
                                    AA_HBondable == TRUE & AA_A_HBondable == FALSE ~ "Added"),
         .before = "Site",
         .keep = "unused")  %>%
  arrange(BondableChange)

writexl::write_xlsx(h_bond, here::here("structure/h_bonds_all.xlsx"))
