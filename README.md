# MSc_Project_Repo
Repository for code and data used for MSc Project 2023/2024: Evolutionary divergence and functional investigation of NAMLAA amidases – bacterial proteins essential to cell division

## Navigation
01_DATA: Store repository for all inputs/outputs for analysis  
03_PROJECT PLANNING: Planning documents and notebook/notes files, and powerpoints for meeting presentations  
04_THESIS: Annotated bibliography and thesis drafts  
05_SCRIPTS: All code and commands used in analysis during this project  
sequence_alignment_files: k-align outputs and cleaned multiple sequence alignment. From previous work, will be deprecated as analysis is updated for 2024  

## Script Index (also listed as Appendix A)

unique_blast_seqs.py: Combine FASTA sequences output from BLASTp searches into a single list and remove duplicate entries based on accession ID.

interpro_download_seqs.py: Download IPR002508 sequences using the InterPro API cropped to the Amidase 3 domain as annotated by InterPro, UniProt accession IDs, taxonomic IDs, and species names.

blast_interpro_collate_annotate.py: Combine sequences from InterPro and BLASTp and remove duplicate sequences using taxonomic ID.

ncbi_taxid_call.py: Annotates BLASTp sequences with the taxonomic ID

call_taxallnomy.py: Annotates sequences with phylum and family based on taxonomic ID.

pre_alignment_filtering.py: Remove sequences with non-specific family names, sequences which failed the API Taxallnomy call, and partial sequences.

gram_stain_updater.R: Filters sequences from the multiple sequence alignment based on information content at each alignment position, and annotates species with predicted gram staining status.

post_alignment_validation.py: Compares structural and sequence alignments to produce True Positive Rate and Positive Predictive Value validation metrics for the multiple sequence alignment.

mobile_helix_conservation.py: Quantifies conservation, identifies insertion regions in the alignment, and produces a binary fingerprint for each sequence based on identified insertion regions.

pca_region_variation.py: Clusters sequences based on binary fingerprint and calculates a 3-component PCA based on fingerprint clusters and scree plot.

domain_structure_prep.py: Downloads all accessions for sequences (NCBI or UniProt).

mobile_helix_conservation_validated_ins.py: Repeats the fingerprinting in mobile_helix_conservation.py with selected reviewed insertion regions.

chimera_code_Amidase 3_structures.txt: Crops experimentally-derived AmiC PDB structures to the Amidase 3 domain.

