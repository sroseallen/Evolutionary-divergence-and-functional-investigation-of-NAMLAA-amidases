if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")

library(tidyr)
library(dplyr)
library(stringr)
library(msa)
library(Biostrings)
library(AMR)

setwd("D:/Documents/MSc_Project/MSc_Project/01_DATA/Amidase_3/04_Multiple_Alignment")

# cleaning of alignment
alignment_2 <- readAAStringSet("alignment_2.fa")
sequence_names <- names(alignment_2)
sequences <- as.character(alignment_2)
alignment_df <- data.frame(sequence_name = sequence_names, sequence = sequences)

alignment_df %>%
  separate_wider_position(col=sequence, c(char=rep(1,1677))) -> alignment_df

# identify gaps inserted by X% of sequences 
## threshold for 0.1% of sequences: a non-gap must be present in at least 40958/1000 sequences (=41 sequences)
## threshold for 1% of sequences: a non-gap must be present in at least 40958/100 sequences (=410 sequences)
gap_threshold = 410

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
  unite(seq, c(2:352), sep="") -> align_fasta

# add AMR annotations to sequence names
cleaner_func <- function(df) {
  df %>%
    mutate (sequence_name := str_replace(sequence_name, ",", "")) %>%
    mutate (sequence_name := str_replace(sequence_name, "\\[", "")) %>%
    mutate (sequence_name := str_replace(sequence_name, "\\]", "")) %>%
    mutate (sequence_name := str_replace(sequence_name, "\\/[[:graph:]][[:graph:]][[:graph:]][[:graph:]][[:graph:]]", "")) %>%
    mutate (sequence_name := str_replace(sequence_name, "\\/[[:graph:]][[:graph:]][[:graph:]][[:graph:]]", "")) %>%
    mutate (sequence_name := str_replace(sequence_name, "\\/[[:graph:]][[:graph:]][[:graph:]]", ""))
}

align_cleaned <- cleaner_func(align_fasta)

# use AMR to add extra information to each species name
AMR_caller <- function (df) {
  df %>%
    mutate(amr_phylum        := mo_phylum       (df$genus_species, keep_synonyms = getOption("AMR_keep_synonyms", FALSE))) %>%
    mutate(amr_fullname      := mo_fullname     (df$genus_species, keep_synonyms = getOption("AMR_keep_synonyms", FALSE))) %>%
    #   mutate(amr_family        := mo_family       (df$genus_species, keep_synonyms = getOption("AMR_keep_synonyms", FALSE))) %>%
    mutate(amr_gram_status   := mo_gramstain    (df$genus_species, keep_synonyms = getOption("AMR_keep_synonyms", FALSE)))
}

align_cleaned %>%
  separate_wider_delim(sequence_name, delim=":", names=c("taxid", "family", "genus_species"), too_many="drop") %>%
  AMR_caller() %>%
  unite(sequence_name, c("taxid", "family", "genus_species", "amr_phylum", "amr_fullname", "amr_gram_status"), sep=":") -> align_fasta_annotated

write.csv(align_fasta_annotated, "temp.csv")
# manual review and update of unknown species
align_fasta_annotated <- read.csv("temp.csv")

# export alignment
AAStringSet(align_fasta_annotated$seq) -> aligned_sequences
align_fasta_annotated$sequence_name -> names(aligned_sequences)
writeXStringSet(aligned_sequences, file="alignment_3_thresh1.0.fa", format="fasta")

# convert cleaned alignment to text file for reference
align_fasta_annotated %>%
  mutate(seq=gsub("-", "", seq)) %>%
  mutate(sequence_name=paste0(">",sequence_name)) -> align_pretxt

txt_file <- file("thresh1.0_seqs.txt", "w")

for (i in 1:nrow(align_pretxt)) {
  writeLines(paste(align_pretxt[i, "sequence_name"]), txt_file)
  writeLines(paste(align_pretxt[i, "seq"]), txt_file)
}

close(txt_file)