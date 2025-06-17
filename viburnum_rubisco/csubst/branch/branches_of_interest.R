library(here)
library(ape)
library(dplyr)
library(ggtree)

print("Reading csubst output")
data <- readr::read_tsv(here("csubst_cb_2.tsv"))
here::i_am("branches_of_interest.R")

filtered <- data %>% filter(omegaCany2spe >= 1 & OCNany2spe >= 1)

branches_of_interest <- filtered %>%
  dplyr::select(branch_id_1, branch_id_2) %>%
  dplyr::add_row(branch_id_1 = 282, branch_id_2 = 293) # nodes before root node
print("Saving branches of interest (ω_C > 1)")
readr::write_csv(branches_of_interest, 
                 file = here("branches_of_interest.txt"),
                 col_names = FALSE)

