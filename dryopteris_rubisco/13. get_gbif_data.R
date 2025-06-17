library(dplyr)
library(purrr)
library(readr)
library(rgbif)
library(CoordinateCleaner)

here::i_am("13. get_gbif_data.R")

# List of species from GenBank download
seq_species <- read.csv(here::here("metadata.csv")) %>% # file with accessions
  mutate(species = stringr::str_remove_all(species, "\\s[:alnum:][:digit:]+$")) %>%
  dplyr::pull(species)

# Find taxon keys for each species in the list
raw_taxon_keys <- map_dfr(seq_species, ~ name_backbone(name = ., verbose = TRUE)) %>% # name_backbone on each species and bind the rows
  filter(!is.na(species)) %>%
  filter(confidence >= 90) %>%
  distinct()

# Split into exact and fuzzy matches, and manually check the fuzzy matches
exact_taxon_keys <- raw_taxon_keys %>%
  filter(matchType == "EXACT")
mismatches <- raw_taxon_keys %>% 
  filter(canonicalName != verbatim_name) %>% 
  select(scientificName, verbatim_name)

FUZZY_SPECIES_KEEP <- c( # Manually check names for reasonable matches, e.g., off by one or two letters
  "Dryopteris affinis subsp. borreri (Newman) Fraser-Jenk.",
  "Dryopteris affinis var. borreri",
  "Dryopteris barbigera (T.Moore ex Hook.) Kuntze",
  "Dryopteris championiae (Benth.) C.Chr.",
  "Dryopteris championiae (Benth.) C.Chr. ex Ching",
  "Dryopteris christenseniae (Ching) Li Bing Zhang",
  "Dryopteris christenseniae (Ching) LiBingZhang",
  "Dryopteris coreanomontana Nakai",
  "Dryopteris fuscoatra var. lamoureuxii Fraser-Jenk.",
  "Dryopteris glabra var. glabra",
  "Dryopteris goldieana (Hook. ex Goldie) A.Gray",
  "Dryopteris insularis subsp. insularis",
  "Dryopteris koidzumiana Tagawa",
  "Dryopteris lunensis (Christ) C.Chr.",
  "Dryopteris simasakii var. simasakii",
  "Dryopteris villarsii Woynar ex Schinz & Thell."
)
fuzzy_taxon_keys <- raw_taxon_keys %>%
  filter(matchType == "FUZZY") %>%
  filter(scientificName %in% FUZZY_SPECIES_KEEP)

# Save taxon keys and species keys (subset of taxon keys)
if (dir.exists("gbif_data") == FALSE) {
  dir.create("gbif_data")
}
taxon_keys <- bind_rows(exact_taxon_keys, fuzzy_taxon_keys) # merge data frames with exact and fuzzy matches
write_csv(taxon_keys, here::here("gbif_data/taxon_keys.csv"))
species_keys <- taxon_keys %>%
  pull(speciesKey)
write_csv(tibble(species_keys), here::here("gbif_data/species_keys.csv"))

# Download occurrences
gbif_download <- occ_download(
  type = "and", 
  pred_in("taxonKey", species_keys),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred_not(pred_in("issue",c("CONTINENT_INVALID",
                             "COUNTRY_COORDINATE_MISMATCH",
                             "COUNTRY_INVALID",
                             "COUNTRY_MISMATCH",
                             "PRESUMED_NEGATED_LATITUDE",
                             "PRESUMED_NEGATED_LONGITUDE",
                             "PRESUMED_SWAPPED_COORDINATE",
                             "TAXON_MATCH_FUZZY",
                             "COORDINATE_PRECISION_INVALID"))),
  pred_gt("decimalLatitude", -60),
  pred_lte("decimalLatitude",90),
  pred_gte("decimalLongitude",-180),
  pred_lte("decimalLongitude",180),
  pred("occurrenceStatus", "PRESENT"),
  format = "SIMPLE_CSV")
occ_download_wait(gbif_download, curlopts=list(http_version=2)) # wait for download to finish
gbif_download <- occ_download_get(gbif_download, path = getwd(), overwrite = TRUE) # save download
gbif_citation <- gbif_citation(gbif_download)
writeLines(gbif_citation$download, "gbif_data/gbif_citation.txt")

gbif_data <- occ_download_import(gbif_download) # load occurrence data

# Write GBIF data into csv
write_csv(gbif_data, here::here("gbif_data/raw.csv"))

gbif_data <- read_csv(here::here("gbif_data/raw.csv"))
taxon_keys <- read_csv(here::here("gbif_data/taxon_keys.csv"))

# Run tests to clean the GBIF data
clean_gbif_data <- gbif_data %>%
  cc_dupl(., lon = "decimalLongitude", lat = "decimalLatitude") %>% # excludes duplicates
  cc_zero(., lon = "decimalLongitude", lat = "decimalLatitude") %>% # excludes zero coordinates
  cc_equ(., lon = "decimalLongitude", lat = "decimalLatitude") %>% # excludes points with identical lat/lon
  cc_sea(., lon = "decimalLongitude", lat = "decimalLatitude", ref = buffland, scale = 10, value = "clean", verbose = TRUE) %>% #excludes oceanic points
  cc_inst(., lon = "decimalLongitude", lat = "decimalLatitude", buffer = 500) %>% # excludes points near biodiversity institutions
  cc_cap(., lon = "decimalLongitude", lat = "decimalLatitude", buffer = 2500) %>% # excludes points near country capitals
  cc_cen(., lon = "decimalLongitude", lat = "decimalLatitude", buffer = 500) %>% # excludes points near country centroids
#  cc_outl(., lon = "decimalLongitude", lat = "decimalLatitude", method = "quantile", mltpl = 5) # excludes geographic outliers
  right_join(x = ., y = dplyr::select(taxon_keys, speciesKey, verbatim_name), by = "speciesKey") # add verbatim_name (the name we searched by) back in
write_csv(clean_gbif_data, "gbif_data/clean.csv")

# Review # of records per species
record_count <- clean_gbif_data %>%
  group_by(verbatim_name) %>%
  distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) %>%
  summarize(RecordCount = n())
write_csv(record_count, "gbif_data/species_sample_counts.csv")
