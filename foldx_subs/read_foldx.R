library(dplyr)
here::i_am("foldx_subs/read_foldx.R")

mutations_list <- readr::read_lines(here::here("foldx_subs/individual_list_mutations.txt")) %>%
  stringr::word(1, sep = ",") %>%
  stringr::str_sub_all(c(1, 3), c(1, -1))

mutations <- purrr::map_chr(mutations_list, ~ paste0(.[1], .[2]))

mutations_seqnum <- tibble::tibble(SeqNum = 1:60, Mutation = mutations)

file_list <- list.files(here::here("foldx_subs/runs"), pattern = "Dif_run")

runs <- readr::read_tsv(here::here("foldx_subs/runs/Dif_run_1_8ruc-assembly-clean_Repair.fxout"), skip = 8) %>%
  mutate(Run = 1 %>% as.numeric(),
         # SeqNum = row_number(),
         SeqNum = rep(1:(n()/3), each = 3),
         .before = "Pdb")

run_num_seq_num <- c("1" = 0,
                     "2" = 8,
                     "3" = 16,
                     "4" = 24,
                     "5" = 32,
                     "6" = 39,
                     "7" = 46,
                     "8" = 53)

for (file in file_list) {
  run_num <- file %>% stringr::str_extract("(?<=run_)[:digit:]+") %>% as.character
  print(paste("Run", run_num))
  if (run_num == "1") {next}
  
  print(paste("Starting at sequence", run_num_seq_num[run_num]))
  run_data <- readr::read_tsv(here::here("foldx_subs", "runs", file), skip = 8) %>%
    filter(stringr::str_detect(Pdb, "^8ruc")) %>%
    mutate(Run = run_num  %>% as.numeric(),
           # SeqNum = row_number() + run_num_seq_num[run_num],
           SeqNum = rep(1:(n()/3), each = 3) + run_num_seq_num[run_num],
           .before = "Pdb")
  prev_seq_num <- run_data %>% pull(SeqNum) %>% max(., na.rm = TRUE)
  runs <- bind_rows(runs, run_data)
}

runs_summary <- runs %>%
  group_by(Run, SeqNum) %>%
  summarise(across(.cols = !c(Pdb),
                   .fns = mean))

runs_sd <- runs %>%
  group_by(Run, SeqNum) %>%
  summarise(SD = sd(`total energy`))

ddG <- runs_summary %>%
  ungroup() %>%
  rename_with(~stringr::str_remove_all(., " ")) %>% 
  mutate(Genus = case_when(SeqNum %>% between(1, 26) ~ "Viburnum",
                           SeqNum %>% between(27, 38) ~ "Dryopteris",
                           SeqNum %>% between(39, 43) ~ "Limonium",
                           SeqNum %>% between(43, 60) ~ "Pinus")) %>%
  full_join(mutations_seqnum, .) %>%
  full_join(runs_sd, .) %>%
  mutate(across(.cols = !c("Genus", "Mutation"),
                .fns = ~ plyr::round_any(., 0.01))) %>%
  mutate(HBond = BackboneHbond + SidechainHbond,
         Clashes = VanderWaalsclashes + torsionalclash + backboneclash,
         .after = totalenergy) %>%
  select(Genus, Mutation, totalenergy, SD, HBond, BackboneHbond, SidechainHbond, Clashes, VanderWaalsclashes, torsionalclash, backboneclash, VanderWaals, Electrostatics, SolvationPolar, SolvationHydrophobic)
  
# Derived AA == spinach AA
##Viburnum
S76N <- ddG %>%
  filter(Genus == "Viburnum", Mutation == "N76S") %>%
  mutate(Mutation = "S76N",
         across(.cols = !c(Genus, Mutation, SD),
                .fns = ~ . * -1))
I225L <- ddG %>%
  filter(Genus == "Viburnum", Mutation == "L225I") %>%
  mutate(Mutation = "I225L",
         across(.cols = !c(Genus, Mutation, SD),
                .fns = ~ . * -1))
E249D <- ddG %>%
  filter(Genus == "Viburnum", Mutation == "D249E") %>%
  mutate(Mutation = "E249D",
         across(.cols = !c(Genus, Mutation, SD),
                .fns = ~ . * -1))
A328S <- ddG %>%
  filter(Genus == "Viburnum", Mutation == "S328A", row_number() == 18) %>%
  mutate(Mutation = "A328S",
         across(.cols = !c(Genus, Mutation, SD),
                .fns = ~ . * -1))
E340D <- ddG %>%# forgot to add this for viburnum, but shouldn't matter if I use limonium one
  filter(Mutation == "D340E") %>%
  mutate(Mutation = "E340D",
         Genus = "Viburnum",
         across(.cols = !c(Genus, Mutation, SD),
                .fns = ~ . * -1))
##Dryopteris
L116M <- ddG %>%
  filter(Genus == "Dryopteris", Mutation == "M116L") %>%
  mutate(Mutation = "L116M",
         across(.cols = !c(Genus, Mutation, SD),
                .fns = ~ . * -1))
A262V <- ddG %>%
  filter(Genus == "Dryopteris", Mutation == "V262A") %>%
  mutate(Mutation = "A262V",
         across(.cols = !c(Genus, Mutation, SD),
                .fns = ~ . * -1))
F282H <- ddG %>%
  filter(Genus == "Dryopteris", Mutation == "H282F") %>%
  mutate(Mutation = "F282H",
         across(.cols = !c(Genus, Mutation, SD),
                .fns = ~ . * -1))
##Limonium
T43S <- ddG %>%
  filter(Genus == "Limonium", Mutation == "S43T") %>%
  mutate(Mutation = "T43S",
         across(.cols = !c(Genus, Mutation, SD),
                .fns = ~ . * -1))
I309M <- ddG %>%
  filter(Genus == "Limonium", Mutation == "M309I") %>%
  mutate(Mutation = "I309M",
         across(.cols = !c(Genus, Mutation, SD),
                .fns = ~ . * -1))
##Pinus
S95N <- ddG %>%
  filter(Genus == "Pinus", Mutation == "N95S") %>%
  mutate(Mutation = "S95N",
         across(.cols = !c(Genus, Mutation, SD),
                .fns = ~ . * -1))
N226Y <- ddG %>%
  filter(Genus == "Pinus", Mutation == "Y226N") %>%
  mutate(Mutation = "N226Y",
         across(.cols = !c(Genus, Mutation, SD),
                .fns = ~ . * -1))
V471A <- ddG %>%
  filter(Genus == "Pinus", Mutation == "A471V") %>%
  mutate(Mutation = "V471A",
         across(.cols = !c(Genus, Mutation, SD),
                .fns = ~ . * -1))
derived_spinach <- bind_rows(S76N, I225L, E249D, A328S, E340D, L116M, A262V, F282H, T43S, I309M, S95N, N226Y, V471A)

# Neither AA == spinach AA
## Viburnum
S95T_v <- ddG %>%
  filter(Genus == "Viburnum", Mutation %in% c("N95S", "N95T")) %>%
  mutate(SD = sum(SD),
         AncDer = ifelse(stringr::str_sub(Mutation, -1, -1) == "T", "derived", "ancestral")) %>% ## DERIVED AA
  tidyr::pivot_wider(names_from = "AncDer",
                     values_from = !c("Genus", "SD")) %>%
  mutate(across(.cols = contains("_derived") & !contains(c("Mutation", "SD", "AncDer")),
                # cur_column gets name of current column, then we replace the _derived suffix. pick uses this string to get the column
                # but pick makes the col a df[1x1], so pluck unnests it
                .fns = ~ . - pick(cur_column() %>% stringr::str_replace("_derived", "_ancestral")) %>% purrr::pluck(1), 
                .names = "diff_{.col}")) %>%
  select(Genus, SD, contains("diff")) %>%
  rename_with(.cols = contains("diff"),
              .fn = ~ stringr::str_extract(., "(?<=\\_)[:graph:]+(?=\\_)")) %>%
  rename_with(~ stringr::str_remove(., "\\$.*")) %>%
  mutate(Mutation = "S95T", .before = "Genus")
A328G_v <- ddG %>%
  filter(Genus == "Viburnum", Mutation %in% c("S328A", "S328G"), row_number() %in% c(19, 20)) %>%
  mutate(SD = sum(SD),
         AncDer = ifelse(stringr::str_sub(Mutation, -1, -1) == "G", "derived", "ancestral")) %>% ## DERIVED AA
  tidyr::pivot_wider(names_from = "AncDer",
                     values_from = !c("Genus", "SD")) %>%
  mutate(across(.cols = contains("_derived") & !contains(c("Mutation", "SD", "AncDer")),
                # cur_column gets name of current column, then we replace the _derived suffix. pick uses this string to get the column
                # but pick makes the col a df[1x1], so pluck unnests it
                .fns = ~ . - pick(cur_column() %>% stringr::str_replace("_derived", "_ancestral")) %>% purrr::pluck(1), 
                .names = "diff_{.col}")) %>%
  select(Genus, SD, contains("diff")) %>%
  rename_with(.cols = contains("diff"),
              .fn = ~ stringr::str_extract(., "(?<=\\_)[:graph:]+(?=\\_)")) %>%
  rename_with(~ stringr::str_remove(., "\\$.*")) %>%
  mutate(Mutation = "A328G", .before = "Genus")
C449A_v <- ddG %>% # forgot to add this for viburnum, but shouldn't matter if I use pinus one
  filter(Genus == "Pinus", Mutation %in% c("T449C", "T449A")) %>%
  mutate(SD = sum(SD),
         AncDer = ifelse(stringr::str_sub(Mutation, -1, -1) == "A", "derived", "ancestral")) %>% ## DERIVED AA
  tidyr::pivot_wider(names_from = "AncDer",
                     values_from = !c("Genus", "SD")) %>%
  mutate(across(.cols = contains("_derived") & !contains(c("Mutation", "SD", "AncDer")),
                # cur_column gets name of current column, then we replace the _derived suffix. pick uses this string to get the column
                # but pick makes the col a df[1x1], so pluck unnests it
                .fns = ~ . - pick(cur_column() %>% stringr::str_replace("_derived", "_ancestral")) %>% purrr::pluck(1), 
                .names = "diff_{.col}")) %>%
  select(Genus, SD, contains("diff")) %>%
  rename_with(.cols = contains("diff"),
              .fn = ~ stringr::str_extract(., "(?<=\\_)[:graph:]+(?=\\_)")) %>%
  rename_with(~ stringr::str_remove(., "\\$.*")) %>%
  mutate(Mutation = "C449A", .before = "Genus") %>%
  mutate(Genus = "Viburnum")
C449S_v <- ddG %>%
  filter(Genus == "Viburnum", Mutation %in% c("T449S", "T449C")) %>%
  mutate(SD = sum(SD),
         AncDer = ifelse(stringr::str_sub(Mutation, -1, -1) == "S", "derived", "ancestral")) %>% ## DERIVED AA
  tidyr::pivot_wider(names_from = "AncDer",
                     values_from = !c("Genus", "SD")) %>%
  mutate(across(.cols = contains("_derived") & !contains(c("Mutation", "SD", "AncDer")),
                # cur_column gets name of current column, then we replace the _derived suffix. pick uses this string to get the column
                # but pick makes the col a df[1x1], so pluck unnests it
                .fns = ~ . - pick(cur_column() %>% stringr::str_replace("_derived", "_ancestral")) %>% purrr::pluck(1), 
                .names = "diff_{.col}")) %>%
  select(Genus, SD, contains("diff")) %>%
  rename_with(.cols = contains("diff"),
              .fn = ~ stringr::str_extract(., "(?<=\\_)[:graph:]+(?=\\_)")) %>%
  rename_with(~ stringr::str_remove(., "\\$.*")) %>%
  mutate(Mutation = "C449S", .before = "Genus")
L475I_v <- ddG %>%
  filter(Genus == "Viburnum", Mutation %in% c("V475L", "V475I")) %>%
  mutate(SD = sum(SD),
         AncDer = ifelse(stringr::str_sub(Mutation, -1, -1) == "I", "derived", "ancestral")) %>% ## DERIVED AA
  tidyr::pivot_wider(names_from = "AncDer",
                     values_from = !c("Genus", "SD")) %>%
  mutate(across(.cols = contains("_derived") & !contains(c("Mutation", "SD", "AncDer")),
                # cur_column gets name of current column, then we replace the _derived suffix. pick uses this string to get the column
                # but pick makes the col a df[1x1], so pluck unnests it
                .fns = ~ . - pick(cur_column() %>% stringr::str_replace("_derived", "_ancestral")) %>% purrr::pluck(1), 
                .names = "diff_{.col}")) %>%
  select(Genus, SD, contains("diff")) %>%
  rename_with(.cols = contains("diff"),
              .fn = ~ stringr::str_extract(., "(?<=\\_)[:graph:]+(?=\\_)")) %>%
  rename_with(~ stringr::str_remove(., "\\$.*")) %>%
  mutate(Mutation = "L475I", .before = "Genus")
## Dryopteris
L116S_d <- ddG %>%
  filter(Genus == "Dryopteris", Mutation %in% c("M116L", "M116S")) %>%
  mutate(SD = sum(SD),
         AncDer = ifelse(stringr::str_sub(Mutation, -1, -1) == "S", "derived", "ancestral")) %>% ## DERIVED AA
  tidyr::pivot_wider(names_from = "AncDer",
                     values_from = !c("Genus", "SD")) %>%
  mutate(across(.cols = contains("_derived") & !contains(c("Mutation", "SD", "AncDer")),
                # cur_column gets name of current column, then we replace the _derived suffix. pick uses this string to get the column
                # but pick makes the col a df[1x1], so pluck unnests it
                .fns = ~ . - pick(cur_column() %>% stringr::str_replace("_derived", "_ancestral")) %>% purrr::pluck(1), 
                .names = "diff_{.col}")) %>%
  select(Genus, SD, contains("diff")) %>%
  rename_with(.cols = contains("diff"),
              .fn = ~ stringr::str_extract(., "(?<=\\_)[:graph:]+(?=\\_)")) %>%
  rename_with(~ stringr::str_remove(., "\\$.*")) %>%
  mutate(Mutation = "L116S", .before = "Genus")
L116P_d <- ddG %>%
  filter(Genus == "Dryopteris", Mutation %in% c("M116L", "M116P")) %>%
  mutate(SD = sum(SD),
         AncDer = ifelse(stringr::str_sub(Mutation, -1, -1) == "P", "derived", "ancestral")) %>% ## DERIVED AA
  tidyr::pivot_wider(names_from = "AncDer",
                     values_from = !c("Genus", "SD")) %>%
  mutate(across(.cols = contains("_derived") & !contains(c("Mutation", "SD", "AncDer")),
                # cur_column gets name of current column, then we replace the _derived suffix. pick uses this string to get the column
                # but pick makes the col a df[1x1], so pluck unnests it
                .fns = ~ . - pick(cur_column() %>% stringr::str_replace("_derived", "_ancestral")) %>% purrr::pluck(1), 
                .names = "diff_{.col}")) %>%
  select(Genus, SD, contains("diff")) %>%
  rename_with(.cols = contains("diff"),
              .fn = ~ stringr::str_extract(., "(?<=\\_)[:graph:]+(?=\\_)")) %>%
  rename_with(~ stringr::str_remove(., "\\$.*")) %>%
  mutate(Mutation = "L116P", .before = "Genus")
L116K_d <- ddG %>%
  filter(Genus == "Dryopteris", Mutation %in% c("M116L", "M116K")) %>%
  mutate(SD = sum(SD),
         AncDer = ifelse(stringr::str_sub(Mutation, -1, -1) == "K", "derived", "ancestral")) %>% ## DERIVED AA
  tidyr::pivot_wider(names_from = "AncDer",
                     values_from = !c("Genus", "SD")) %>%
  mutate(across(.cols = contains("_derived") & !contains(c("Mutation", "SD", "AncDer")),
                # cur_column gets name of current column, then we replace the _derived suffix. pick uses this string to get the column
                # but pick makes the col a df[1x1], so pluck unnests it
                .fns = ~ . - pick(cur_column() %>% stringr::str_replace("_derived", "_ancestral")) %>% purrr::pluck(1), 
                .names = "diff_{.col}")) %>%
  select(Genus, SD, contains("diff")) %>%
  rename_with(.cols = contains("diff"),
              .fn = ~ stringr::str_extract(., "(?<=\\_)[:graph:]+(?=\\_)")) %>%
  rename_with(~ stringr::str_remove(., "\\$.*")) %>%
  mutate(Mutation = "L116K", .before = "Genus")
F282Y_d <- ddG %>%
  filter(Genus == "Dryopteris", Mutation %in% c("H282F", "H282Y")) %>%
  mutate(SD = sum(SD),
         AncDer = ifelse(stringr::str_sub(Mutation, -1, -1) == "Y", "derived", "ancestral")) %>% ## DERIVED AA
  tidyr::pivot_wider(names_from = "AncDer",
                     values_from = !c("Genus", "SD")) %>%
  mutate(across(.cols = contains("_derived") & !contains(c("Mutation", "SD", "AncDer")),
                # cur_column gets name of current column, then we replace the _derived suffix. pick uses this string to get the column
                # but pick makes the col a df[1x1], so pluck unnests it
                .fns = ~ . - pick(cur_column() %>% stringr::str_replace("_derived", "_ancestral")) %>% purrr::pluck(1), 
                .names = "diff_{.col}")) %>%
  select(Genus, SD, contains("diff")) %>%
  rename_with(.cols = contains("diff"),
              .fn = ~ stringr::str_extract(., "(?<=\\_)[:graph:]+(?=\\_)")) %>%
  rename_with(~ stringr::str_remove(., "\\$.*")) %>%
  mutate(Mutation = "F282Y", .before = "Genus")
F282L_d <- ddG %>%
  filter(Genus == "Dryopteris", Mutation %in% c("H282F", "H282L")) %>%
  mutate(SD = sum(SD),
         AncDer = ifelse(stringr::str_sub(Mutation, -1, -1) == "L", "derived", "ancestral")) %>% ## DERIVED AA
  tidyr::pivot_wider(names_from = "AncDer",
                     values_from = !c("Genus", "SD")) %>%
  mutate(across(.cols = contains("_derived") & !contains(c("Mutation", "SD", "AncDer")),
                # cur_column gets name of current column, then we replace the _derived suffix. pick uses this string to get the column
                # but pick makes the col a df[1x1], so pluck unnests it
                .fns = ~ . - pick(cur_column() %>% stringr::str_replace("_derived", "_ancestral")) %>% purrr::pluck(1), 
                .names = "diff_{.col}")) %>%
  select(Genus, SD, contains("diff")) %>%
  rename_with(.cols = contains("diff"),
              .fn = ~ stringr::str_extract(., "(?<=\\_)[:graph:]+(?=\\_)")) %>%
  rename_with(~ stringr::str_remove(., "\\$.*")) %>%
  mutate(Mutation = "F282L", .before = "Genus")
F282I_d <- ddG %>%
  filter(Genus == "Dryopteris", Mutation %in% c("H282F", "H282I")) %>%
  mutate(SD = sum(SD),
         AncDer = ifelse(stringr::str_sub(Mutation, -1, -1) == "I", "derived", "ancestral")) %>% ## DERIVED AA
  tidyr::pivot_wider(names_from = "AncDer",
                     values_from = !c("Genus", "SD")) %>%
  mutate(across(.cols = contains("_derived") & !contains(c("Mutation", "SD", "AncDer")),
                # cur_column gets name of current column, then we replace the _derived suffix. pick uses this string to get the column
                # but pick makes the col a df[1x1], so pluck unnests it
                .fns = ~ . - pick(cur_column() %>% stringr::str_replace("_derived", "_ancestral")) %>% purrr::pluck(1), 
                .names = "diff_{.col}")) %>%
  select(Genus, SD, contains("diff")) %>%
  rename_with(.cols = contains("diff"),
              .fn = ~ stringr::str_extract(., "(?<=\\_)[:graph:]+(?=\\_)")) %>%
  rename_with(~ stringr::str_remove(., "\\$.*")) %>%
  mutate(Mutation = "F282I", .before = "Genus")
## Pinus
S95T_p <- ddG %>%
  filter(Genus == "Pinus", Mutation %in% c("N95S", "N95T")) %>%
  mutate(SD = sum(SD),
         AncDer = ifelse(stringr::str_sub(Mutation, -1, -1) == "T", "derived", "ancestral")) %>% ## DERIVED AA
  tidyr::pivot_wider(names_from = "AncDer",
                     values_from = !c("Genus", "SD")) %>%
  mutate(across(.cols = contains("_derived") & !contains(c("Mutation", "SD", "AncDer")),
                # cur_column gets name of current column, then we replace the _derived suffix. pick uses this string to get the column
                # but pick makes the col a df[1x1], so pluck unnests it
                .fns = ~ . - pick(cur_column() %>% stringr::str_replace("_derived", "_ancestral")) %>% purrr::pluck(1), 
                .names = "diff_{.col}")) %>%
  select(Genus, SD, contains("diff")) %>%
  rename_with(.cols = contains("diff"),
              .fn = ~ stringr::str_extract(., "(?<=\\_)[:graph:]+(?=\\_)")) %>%
  rename_with(~ stringr::str_remove(., "\\$.*")) %>%
  mutate(Mutation = "S95T", .before = "Genus")
S95A_p <- ddG %>%
  filter(Genus == "Pinus", Mutation %in% c("N95S", "N95A")) %>%
  mutate(SD = sum(SD),
         AncDer = ifelse(stringr::str_sub(Mutation, -1, -1) == "A", "derived", "ancestral")) %>% ## DERIVED AA
  tidyr::pivot_wider(names_from = "AncDer",
                     values_from = !c("Genus", "SD")) %>%
  mutate(across(.cols = contains("_derived") & !contains(c("Mutation", "SD", "AncDer")),
                # cur_column gets name of current column, then we replace the _derived suffix. pick uses this string to get the column
                # but pick makes the col a df[1x1], so pluck unnests it
                .fns = ~ . - pick(cur_column() %>% stringr::str_replace("_derived", "_ancestral")) %>% purrr::pluck(1), 
                .names = "diff_{.col}")) %>%
  select(Genus, SD, contains("diff")) %>%
  rename_with(.cols = contains("diff"),
              .fn = ~ stringr::str_extract(., "(?<=\\_)[:graph:]+(?=\\_)")) %>%
  rename_with(~ stringr::str_remove(., "\\$.*")) %>%
  mutate(Mutation = "S95A", .before = "Genus")
C449A_p <- ddG %>%
  filter(Genus == "Pinus", Mutation %in% c("T449C", "T449A")) %>%
  mutate(SD = sum(SD),
         AncDer = ifelse(stringr::str_sub(Mutation, -1, -1) == "A", "derived", "ancestral")) %>% ## DERIVED AA
  tidyr::pivot_wider(names_from = "AncDer",
                     values_from = !c("Genus", "SD")) %>%
  mutate(across(.cols = contains("_derived") & !contains(c("Mutation", "SD", "AncDer")),
                # cur_column gets name of current column, then we replace the _derived suffix. pick uses this string to get the column
                # but pick makes the col a df[1x1], so pluck unnests it
                .fns = ~ . - pick(cur_column() %>% stringr::str_replace("_derived", "_ancestral")) %>% purrr::pluck(1), 
                .names = "diff_{.col}")) %>%
  select(Genus, SD, contains("diff")) %>%
  rename_with(.cols = contains("diff"),
              .fn = ~ stringr::str_extract(., "(?<=\\_)[:graph:]+(?=\\_)")) %>%
  rename_with(~ stringr::str_remove(., "\\$.*")) %>%
  mutate(Mutation = "C449A", .before = "Genus")
C449S_p <- ddG %>%
  filter(Genus == "Pinus", Mutation %in% c("T449C", "T449S")) %>%
  mutate(SD = sum(SD),
         AncDer = ifelse(stringr::str_sub(Mutation, -1, -1) == "S", "derived", "ancestral")) %>% ## DERIVED AA
  tidyr::pivot_wider(names_from = "AncDer",
                     values_from = !c("Genus", "SD")) %>%
  mutate(across(.cols = contains("_derived") & !contains(c("Mutation", "SD", "AncDer")),
                # cur_column gets name of current column, then we replace the _derived suffix. pick uses this string to get the column
                # but pick makes the col a df[1x1], so pluck unnests it
                .fns = ~ . - pick(cur_column() %>% stringr::str_replace("_derived", "_ancestral")) %>% purrr::pluck(1), 
                .names = "diff_{.col}")) %>%
  select(Genus, SD, contains("diff")) %>%
  rename_with(.cols = contains("diff"),
              .fn = ~ stringr::str_extract(., "(?<=\\_)[:graph:]+(?=\\_)")) %>%
  rename_with(~ stringr::str_remove(., "\\$.*")) %>%
  mutate(Mutation = "C449S", .before = "Genus")
neither_spinach <- bind_rows(S95T_v, A328G_v, C449A_v, C449S_v, L475I_v, L116S_d, L116P_d, L116K_d, F282Y_d, F282L_d, F282I_d, S95T_p, S95A_p, C449A_p, C449S_p)

# combine all
subs_v <- readr::read_csv(here::here("analyze/aa_at_sites_ancestral.csv")) %>%
  mutate(Genus = "Viburnum") 
subs_d <- readr::read_csv(here::here("dryopteris_rubisco/aa_at_sites_ancestral.csv")) %>%
  mutate(Genus = "Dryopteris") 
subs_l <- readr::read_csv(here::here("limonium_rubisco/aa_at_sites_ancestral.csv")) %>%
  mutate(Genus = "Limonium") 
subs_p <- readr::read_csv(here::here("pine_rubisco/aa_at_sites_ancestral.csv")) %>%
  mutate(Genus = "Pinus")
subs <- bind_rows(subs_v, subs_d, subs_l, subs_p) %>%
  rowwise() %>%
  filter(!AA %in% c(AA_Ancestral, "-")) %>%
  distinct(Genus, Site, AA_Ancestral, Site, AA) %>%
  mutate(Substitution = paste0(AA_Ancestral, Site, AA)) %>%
  select(!contains("AA"))

ddG_all <- bind_rows(ddG, derived_spinach, neither_spinach) %>%
  rename(Substitution = Mutation) %>%
  right_join(subs) %>%
  relocate(Site, .before = "Substitution")

readr::write_csv(ddG_all, here::here("foldx_subs/ddG_values.csv"))                
