library(dplyr)

here::i_am("12. merge_site_data.R")

# directory list have fixed (must match fully) pattern of "csubst_site.branch_id"

dir_list <- dir(here::here("site"), pattern = "^csubst_site.branch_id")


# function that splits directory name to get the branch numbers, then reads each tsv in each directory. puts in the branch numbers into the data frame before returning it

read_csubst_site <- function(dirname) {
  branches <- stringr::str_extract(dirname, "(?<=_id).*")
  branch_1 <- stringr::word(branches, sep = ",", 1)
  branch_2 <- stringr::word(branches, sep = ",", 2)
  print(paste(branch_1, branch_2))
  data <- readr::read_tsv(here::here("site", dirname, "csubst_site.tsv")) %>%
    mutate(branch_id = branches, .before = "codon_site_alignment") %>%
    tidyr::separate(col = branch_id, into = c("branch_id_1", "branch_id_2"), sep = ",") %>%
    rename_with(~ stringr::str_replace(., paste0("_", branch_1, "(?!\\d)"), "_1"), !contains("branch_id")) %>%
    rename_with(~ stringr::str_replace(., paste0("_", branch_2, "(?!\\d)"), "_2"), !contains("branch_id"))
  return(data)
}

test_file <- dir_list[1]
test_file <- "csubst_site.branch_id7,175"
# test_file <- "csubst_site.branch_id10,114"
test_data <- read_csubst_site(test_file)

# run the function across all directories in the site.

merged <- purrr::map_dfr(dir_list, ~ read_csubst_site(.))

# get taxa from each branch

csubst_tree <- ape::read.tree(here::here("branch/csubst_tree.nwk"))
sessa_tree <- ape::read.tree(here::here("dryopteris_csubst.tre"))

branches <- readr::read_tsv(here::here("branch/csubst_b.tsv")) %>%
  rowwise() %>%
  mutate(
    branch_name = ifelse(
      stringr::str_detect(branch_name, "^Node"), #if contains Node
      purrr::keep(csubst_tree$node.label, .p = stringr::str_detect, pattern = stringr::fixed(paste(branch_name, "|", sep =""))), # look in the tree node labels for a partial match
      branch_name 
      )
    )

branches_tips <- branches %>%
  filter(!stringr::str_detect(branch_name, "^Node")) %>%
  mutate(taxa = branch_name) %>% # otherwise just branch_name as taxon name
  dplyr::select(branch_name, branch_id, taxa) %>%
  rowwise() %>%
  mutate(node = branch_name)
branches_nodes <- branches %>%
  filter(stringr::str_detect(branch_name, "^Node")) %>%
  mutate(taxa = list(ape::extract.clade(csubst_tree, branch_name)$tip.label)) %>%
  group_by(branch_name, branch_id) %>%
  summarise(taxa = purrr::map(taxa, ~ stringr::word(., 1, sep = stringr::fixed("|")))) %>%
  mutate(taxa = purrr::list_merge(taxa)) %>%
  rowwise() %>%
  mutate(node = list(ape::getMRCA(sessa_tree, taxa))) %>% # get node numbers from Jin 2021 PNAS tree
  tidyr::unnest(node)

branches_all <- rbind(branches_tips, branches_nodes) %>%
  dplyr::select(!taxa)

readr::write_csv(branches_all, here::here("branch/branch_id_nodes.csv"))

# add node numbers to merged data

node1 <- branches_all %>% 
  rename(branch_id_1 = branch_id, node1 = node) %>%
  mutate(branch_id_1 = as.character(branch_id_1))
node2 <- branches_all %>%
  rename(branch_id_2 = branch_id, node2 = node) %>%
  mutate(branch_id_2 = as.character(branch_id_2))

merged_nodes <- full_join(merged, node1, by = "branch_id_1") %>%
  full_join(node2, by = "branch_id_2") %>%
  relocate(c(node1, node2), .after = "branch_id_2")
# save

readr::write_csv(merged_nodes, here::here("site/csubst_site_merged.csv"))

