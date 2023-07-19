library(AMR)
options(AMR_keep_synonyms = TRUE)

library(tidyverse)

setwd("MSc_Project/MSc_Project/")

# read in file with bacterial names
zero <- read_csv("0_percent.csv")
fifty <- read_csv("50_percent.csv")
seventy <- read_csv("70_percent.csv")
ninety <- read_csv("90_percent.csv")

# clean up genus/species column for running in AMR package
cleaner_func <- function(df) {
  df %>% 
    select (-1) %>%
    mutate (genus_species := str_replace(genus_species, ",", "")) %>%
    mutate (genus_species := str_replace(genus_species, "\\[", "")) %>%
    mutate (genus_species := str_replace(genus_species, "\\]", "")) %>%
    mutate (genus_species := str_replace(genus_species, "\\/[[:graph:]][[:graph:]][[:graph:]][[:graph:]][[:graph:]]", "")) %>%
    mutate (genus_species := str_replace(genus_species, "\\/[[:graph:]][[:graph:]][[:graph:]][[:graph:]]", "")) %>%
    mutate (genus_species := str_replace(genus_species, "\\/[[:graph:]][[:graph:]][[:graph:]]", ""))
}

cleaner_func(zero) -> zero
cleaner_func(fifty) -> fifty
cleaner_func(seventy) -> seventy
cleaner_func(ninety) -> ninety


# AMR for information about each bacterial species
AMR_caller <- function (df) {
  df %>%
    mutate(amr_fullname      := mo_fullname     (df$genus_species)) %>%
    mutate(amr_phylum        := mo_phylum       (df$genus_species)) %>%
    mutate(amr_family        := mo_family       (df$genus_species)) %>%
    mutate(amr_pathogenicity := mo_pathogenicity(df$genus_species)) %>%
    mutate(amr_gram_status   := mo_gramstain    (df$genus_species))
}

AMR_caller(zero) -> zero
AMR_caller(fifty) -> fifty
AMR_caller(seventy) -> seventy
AMR_caller(ninety) -> ninety

# save results to csv
write_csv(zero, "0_percent_info.csv")
write_csv(fifty, "50_percent_info.csv")
write_csv(seventy, "70_percent_info.csv")
write_csv(ninety, "90_percent_info.csv")
