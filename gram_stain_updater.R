if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("odseq")

library(AMR)
library(tidyverse)
library(msa)
library(odseq)

setwd("MSc_Project/MSc_Project/")

# Step 1: KALIGN + OD-SEQ

# Step 2: Cleaning of alignment, phase 2:

read_csv("MSA_cleaning/affine_core_1_stdev_clean.csv", col_names=FALSE) -> alignment

alignment %>%
  separate(col = X2, into = paste0("char_", 1:1395), sep = "") %>%
  select (-2) -> align_df

# identify gaps inserted by X% of sequences (29 sequences max. involved in an insertion position for 0.1%)

gap_threshold = 286

align_df %>% top_n(0) -> align_removed
(colnames(align_df))[-1]-> align_cols

for (col_name in align_cols) {
  align_df %>% filter(across(all_of(col_name), ~ . != "-")) -> temp
  if (nrow(temp) <= gap_threshold) {
    bind_rows(temp, align_removed) -> align_removed
    align_df %>% filter (across(all_of(col_name), ~ . == "-")) -> align_df
  }
}

# remove columns of only gaps
align_df %>%
  select(where(~ any(. != "-"))) -> align_df

# convert cleaned alignment to a fasta file
align_df %>%
  unite(seq, c(2:), sep="", remove=TRUE) -> align_fasta

AAStringSet(align_fasta$seq) -> aligned_sequences
align_fasta$X1 -> names(aligned_sequences)
writeXStringSet(aligned_sequences, file="affine_core_1_stdev_clean_phase2_1.0.fa", format="fasta")

# Step 3: PYTHON. Script 'annotate_seqs_taxid' to get region percent files from the alignment

# Step 4: Read in file with bacterial names
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
    mutate(amr_fullname      := mo_fullname     (df$genus_species, keep_synonyms = getOption("AMR_keep_synonyms", FALSE))) %>%
    mutate(amr_phylum        := mo_phylum       (df$genus_species, keep_synonyms = getOption("AMR_keep_synonyms", FALSE))) %>%
    mutate(amr_family        := mo_family       (df$genus_species, keep_synonyms = getOption("AMR_keep_synonyms", FALSE))) %>%
    mutate(amr_pathogenicity := mo_pathogenicity(df$genus_species, keep_synonyms = getOption("AMR_keep_synonyms", FALSE))) %>%
    mutate(amr_gram_status   := mo_gramstain    (df$genus_species, keep_synonyms = getOption("AMR_keep_synonyms", FALSE)))
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
