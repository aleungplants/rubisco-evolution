library(here)
library(ape)
library(dplyr)
library(ggtree)

here::i_am("10. branches_of_interest.R")

print("Reading csubst output")
data <- readr::read_tsv(here("branch/csubst_cb_2.tsv"))

filtered <- data %>% filter(omegaCany2spe >= 1 & OCNany2spe >= 1)

branches_of_interest <- filtered %>%
  dplyr::select(branch_id_1, branch_id_2)
print("Saving branches of interest (ω_C > 1)")
readr::write_csv(branches_of_interest, 
                 file = here("branch/branches_of_interest.txt"),
                 col_names = FALSE)

