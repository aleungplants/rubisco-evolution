library(dplyr)

downloaded_species <- read.csv("/Users/art/Library/CloudStorage/OneDrive-UniversityofToronto/Viburnum\ rubisco/download/accessions.csv") %>% rename(d_species = species)

landis_species <- read.csv("/Users/art/Library/CloudStorage/OneDrive-UniversityofToronto/Viburnum\ rubisco/download/landis_species.csv", header = FALSE) %>% rename(l_species = V1)

d_missing <- downloaded_species %>% 
  mutate(missing = ifelse(!d_species %in% landis_species$l_species, d_species, NA)) %>%
  select(missing) %>%
  na.omit() %>%
  mutate(source = "download")

l_missing <- landis_species %>% 
  mutate(missing = ifelse(!l_species %in% downloaded_species$d_species, l_species, NA)) %>%
  select(missing) %>%
  na.omit() %>%
  mutate(source = "landis")

missing <- rbind(d_missing, l_missing) %>%
  arrange(missing)

write.csv(missing, here::here("missing_from_download.csv"), row.names = FALSE)
