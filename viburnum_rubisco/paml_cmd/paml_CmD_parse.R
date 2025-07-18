library(dplyr)

here::i_am("paml_cmc/paml_CmC_parse.R")

# PARSE FOR LRT CALCULATIONS AND PARAMETER REPORTING

getAIC <- function(lnL, nParam) {return(-2*(lnL)+2*nParam)}

mlc_CmC <- readr::read_file(here::here("paml_cmc/mlc"))
mlc_M3 <- readr::read_file(here::here("paml/mlc"))

lnL_CmC <- mlc_CmC %>% 
  stringr::str_extract("lnL(.*)") %>%
  stringr::str_extract("[:punct:][:digit:]+[:punct:][:digit:]+") %>%
  as.numeric()
lnL_M3 <- mlc_M3 %>% 
  stringr::str_extract("lnL(.*)") %>%
  stringr::str_extract("[:punct:][:digit:]+[:punct:][:digit:]+") %>%
  as.numeric()

np_CmC <- mlc_CmC %>% stringr::str_extract("(?<=np:)[:digit:]+") %>%
  as.numeric()
np_M3 <- mlc_M3 %>% stringr::str_extract("(?<=np:)[:digit:]+") %>%
  as.numeric()


kappa_CmC <- mlc_CmC %>% 
  stringr::str_extract("kappa(.*)") %>%
  stringr::word(-1)
kappa_M3 <- mlc_M3 %>% 
  stringr::str_extract("kappa(.*)") %>%
  stringr::word(-1)

props_CmC <- mlc_CmC %>% 
  stringr::str_extract("proportion(.*)") %>%
  stringr::str_squish()
props_M3 <- mlc_M3 %>% 
  stringr::str_extract("\np:(.*)") %>%
  stringr::str_squish()

omegas_CmC_bg <- mlc_CmC %>% 
  stringr::str_extract("branch type 0(.*)") %>%
  stringr::str_squish()
omegas_CmC_fg1 <- mlc_CmC %>% 
  stringr::str_extract("branch type 1(.*)") %>%
  stringr::str_squish()
omegas_CmC_fg2 <- mlc_CmC %>% 
  stringr::str_extract("branch type 2(.*)") %>%
  stringr::str_squish()
omegas_CmC_fg3 <- mlc_CmC %>% 
  stringr::str_extract("branch type 3(.*)") %>%
  stringr::str_squish()
omegas_M3 <- mlc_M3 %>% 
  stringr::str_extract("\nw:(.*)") %>%
  stringr::str_squish()

# Model-null model pairs and their comparisons are manually entered
lrt <- tibble(Model = "M3",
              np = np_M3,
              AIC = getAIC(lnL_M3, np_M3),
              lnL = lnL_M3,
              K = kappa_M3,
              w0 = stringr::word(omegas_M3, -3),
              w1 = stringr::word(omegas_M3, -2),
              w2 = stringr::word(omegas_M3, -1),
              p0 = stringr::word(props_M3, -3),
              p1 = stringr::word(props_M3, -2),
              p2 = stringr::word(props_M3, -1),
              ChiSq = NA,
              df = NA,
              p = NA) %>%
  add_row(Model = "BG: warm temperate",
          np = np_CmC,
          AIC = getAIC(lnL_CmC, np_CmC),
          lnL = lnL_CmC,
          K = kappa_CmC,
          w0 = stringr::word(omegas_CmC_bg, -3),
          w1 = stringr::word(omegas_CmC_bg, -2),
          w2 = stringr::word(omegas_CmC_bg, -1),
          p0 = stringr::word(props_CmC, -3),
          p1 = stringr::word(props_CmC, -2),
          p2 = stringr::word(props_CmC, -1),
          ChiSq = 2*(lnL_CmC-lnL_M3),
          df = np_CmC-np_M3,
          p = pchisq(ChiSq, df, lower.tail = FALSE)) %>%
  ## Add paramters
  add_row(Model = "FG: cold temperate",
          w2 = stringr::word(omegas_CmC_fg1, -1)) %>%
  add_row(Model = "FG: tropical",
          w2 = stringr::word(omegas_CmC_fg2, -1)) %>%
  add_row(Model = "FG: cloud",
          w2 = stringr::word(omegas_CmC_fg3, -1)) %>%
  ## Clean up decimals
  mutate(across(.cols = starts_with("w"),
                .fns = ~ ifelse(is.na(.), ., as.numeric(.) %>% plyr::round_any(0.1) %>% as.character()))) %>%
  mutate(across(.cols = c(p0, p1, p2),
                .fns = ~ ifelse(is.na(.), ., (as.numeric(.) * 100) %>% plyr::round_any(0.1) %>% as.character()))) %>%
  mutate(lnL = plyr::round_any(as.numeric(lnL), 0.1) %>% as.character,
         ChiSq = plyr::round_any(as.numeric(ChiSq), 0.1) %>% as.character,
         AIC = plyr::round_any(as.numeric(AIC), 1) %>% as.character,
         K = plyr::round_any(as.numeric(K), 0.01) %>% as.character,
         p = plyr::round_any(as.numeric(p), 0.001) %>% as.character) %>%
  ## Combine proportions with omegas
  mutate(w0 = ifelse(is.na(w0), w0, paste0(w0, " (", p0, "%)")),
         w1 = ifelse(is.na(w0), w1, paste0(w1, " (", p1, "%)")),
         w2 = ifelse(is.na(w0), w2, paste0(w2, " (", p2, "%)")),
         .keep = "unused") %>%
  ## Clean up column names
  rename(`n.p.` = np,
         `κ` = K,
         `χ2` = ChiSq)

writexl::write_xlsx(lrt, here::here("table_PAML_CmC_LRTs.xlsx"))

