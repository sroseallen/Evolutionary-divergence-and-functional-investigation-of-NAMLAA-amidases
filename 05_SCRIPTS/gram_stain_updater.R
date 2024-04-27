if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")

library(tidyr)
library(dplyr)
library(msa)
library(Biostrings)
library(AMR)

setwd("D:/Documents/MSc_Project/MSc_Project/01_DATA/Amidase_3/04_Multiple_Alignment")

# Step 1: KALIGN + OD-SEQ

# Step 2: Cleaning of alignment, phase 2:

alignment_2 <- readAAStringSet("alignment_2.fa")
sequence_names <- names(alignment_2)
sequences <- as.character(alignment_2)
alignment_df <- data.frame(sequence_name = sequence_names, sequence = sequences)

alignment_df %>%
  separate_wider_position(col=sequence, c(char=rep(1,1765))) -> alignment_df

# identify gaps inserted by X% of sequences 
## threshold for 1% of sequences: a non-gap must be present in at least 42660/100 sequences (=427 sequences)

gap_threshold = 427

alignment_df %>% top_n(0) -> align_removed
(colnames(alignment_df))[-1]-> align_cols

for (col_name in align_cols) {
  # identify how many sequences for each column provide data (not a gap)
  alignment_df %>% filter(across(all_of(col_name), ~ . != "-")) -> interim_df
  # if there are fewer of these sequences than the threshold, then that column contains <1% information and all these sequences should be removed.
  if (nrow(interim_df) <= gap_threshold) {
    bind_rows(interim_df, align_removed) -> align_removed
    # keep the remaining sequences, recycle the dataframe to look at the next column
    alignment_df %>% filter (across(all_of(col_name), ~ . == "-")) -> alignment_df
  }
}

## remove columns of only gaps
alignment_df %>%
  select(where(~ any(. != "-"))) -> alignment_keep

# convert cleaned alignment to a fasta file
alignment_keep %>%
  unite(seq, c(2:377), sep="") -> align_fasta

AAStringSet(align_fasta$seq) -> aligned_sequences
align_fasta$sequence_name -> names(aligned_sequences)
writeXStringSet(aligned_sequences, file="alignment_3.fa", format="fasta")

# convert cleaned alignment to text file for re-alignment
align_fasta %>%
  mutate(seq=gsub("-", "", seq)) %>%
  mutate(sequence_name=paste0(">",sequence_name)) -> align_pretxt

txt_file <- file("post_filtering_for_realignment.txt", "w")

for (i in 1:nrow(align_pretxt)) {
  writeLines(paste(align_pretxt[i, "sequence_name"]), txt_file)
  writeLines(paste(align_pretxt[i, "seq"]), txt_file)
}

close(txt_file)

# Step 3: PYTHON. Script 'annotate_seqs_taxid' to get region percent files from the alignment

# Step 4: Read in file with bacterial names
# zero <- read_csv("0_percent.csv")
# fifty <- read_csv("50_percent.csv")
# seventy <- read_csv("70_percent.csv")
# ninety <- read_csv("90_percent.csv")
# 
# # clean up genus/species column for running in AMR package
# cleaner_func <- function(df) {
#   df %>% 
#     select (-1) %>%
#     mutate (genus_species := str_replace(genus_species, ",", "")) %>%
#     mutate (genus_species := str_replace(genus_species, "\\[", "")) %>%
#     mutate (genus_species := str_replace(genus_species, "\\]", "")) %>%
#     mutate (genus_species := str_replace(genus_species, "\\/[[:graph:]][[:graph:]][[:graph:]][[:graph:]][[:graph:]]", "")) %>%
#     mutate (genus_species := str_replace(genus_species, "\\/[[:graph:]][[:graph:]][[:graph:]][[:graph:]]", "")) %>%
#     mutate (genus_species := str_replace(genus_species, "\\/[[:graph:]][[:graph:]][[:graph:]]", ""))
# }
# 
# cleaner_func(zero) -> zero
# cleaner_func(fifty) -> fifty
# cleaner_func(seventy) -> seventy
# cleaner_func(ninety) -> ninety
# 
# # AMR for information about each bacterial species
# AMR_caller <- function (df) {
#   df %>%
#     mutate(amr_fullname      := mo_fullname     (df$genus_species, keep_synonyms = getOption("AMR_keep_synonyms", FALSE))) %>%
#     mutate(amr_phylum        := mo_phylum       (df$genus_species, keep_synonyms = getOption("AMR_keep_synonyms", FALSE))) %>%
#     mutate(amr_family        := mo_family       (df$genus_species, keep_synonyms = getOption("AMR_keep_synonyms", FALSE))) %>%
#     mutate(amr_pathogenicity := mo_pathogenicity(df$genus_species, keep_synonyms = getOption("AMR_keep_synonyms", FALSE))) %>%
#     mutate(amr_gram_status   := mo_gramstain    (df$genus_species, keep_synonyms = getOption("AMR_keep_synonyms", FALSE)))
# }
# 
# AMR_caller(zero) -> zero
# AMR_caller(fifty) -> fifty
# AMR_caller(seventy) -> seventy
# AMR_caller(ninety) -> ninety
# 
# # save results to csv
# write_csv(zero, "0_percent_info.csv")
# write_csv(fifty, "50_percent_info.csv")
# write_csv(seventy, "70_percent_info.csv")
# write_csv(ninety, "90_percent_info.csv")
