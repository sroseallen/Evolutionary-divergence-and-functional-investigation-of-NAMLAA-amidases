from Bio import SeqIO
from Bio import AlignIO
import requests
from tqdm import tqdm
import pandas as pd

# next: write python script to filter out 
#non-specific family names (such as ‘multi-family’)
#non-bacterial Thermococcaceae, 
#and any sequence annotated as being partial

# this script should also return the number of sequences which were not annotated by Taxallnomy, 
#which should be 0 as I have taxids for all of them now?

# N = 56,794
# Removed X because Taxallnomy call failed; N = 56,794 - X

# Removed: Archea, Eukaryota, Virus (non-bacterial kingdom sequences)
## WP_048045911.1:36-222_N-acetylmuramoyl-L-alanine_amidase_[Methanosarcina_mazei] (archea that produces methane, removed)
## Thermococcaceae x29 (archea, extreme thermophiles, removed)
# New N = 

# removed partial sequences
# New N = 



# tag sequences with identity in the helix 5 region
# Based on predicted helix region from 3NE8 the sequence pattern is DAIAKSLAESENKVDLLDG.....DILLDLTRRET. treat these as helical positions, anything that aligns in >90% is high confirmation, >70% as medium confirmation, >50% as low confirmation
# In current alignment: Helix 1: 754, 757, 760, 767, 769, 771, 775, 777, 780, 789, 791, 792, 799, 803, 806, 810, 815, 821, 823. Helix 2: 897, 900, 901, 902, 903, 905, 907, 908, 911, 913, 917

"""
from Bio import AlignIO
alignment = AlignIO.read("./Post-cleaned MFA analyses/affine_core_1_stdev_clean_phase2_1.0.fa", "fasta")

helix_predicted_at = list(range(140, 215)) # length = 75

df = pd.DataFrame({
    "taxid": [],
    "family": [],
    "genus_species":[],
    "helix_5_occupancy_75_max": []
})

for seq in alignment:
    seq_occupancy = 0
    seq_split = seq.description.split(":")
    for pos in helix_predicted_at:
        if seq[pos-1] != '-':
            seq_occupancy += 1
    seq_split.append(seq_occupancy)
    seq_split = pd.DataFrame([seq_split], columns=df.columns)
    df = pd.concat([df, seq_split], ignore_index=True)
    
pd90 = df[df["helix_5_occupancy_75_max"] >= 67]  
pd90.to_csv("90_percent.csv")

pd70 = df[df["helix_5_occupancy_75_max"] >= 52] 
pd70.to_csv("70_percent.csv")

pd50 = df[df["helix_5_occupancy_75_max"] >= 37] 
pd50.to_csv("50_percent.csv")

pd0 = df[df["helix_5_occupancy_75_max"] < 36] 
pd0.to_csv("0_percent.csv")


# % of identicial/consistent class/phylum in the tagged sequences? In the untagged sequences?
## charts made, not really any similarity at family level.
## output taxid, run through taxallnomy, any similarity at hiogher level? 
pd90["taxid"].to_csv("list_90_tax.txt", index=False, header=False, sep=",")
pd70["taxid"].to_csv("list_70_tax.txt", index=False, header=False, sep=",")
pd50["taxid"].to_csv("list_50_tax.txt", index=False, header=False, sep=",")
pd0["taxid"].to_csv("list_0_tax.txt", index=False, header=False, sep=",")

# 90% group contains 15 bacteriophage sequences and 4 gram positive bacterial sequences.
# no mycobacterium phylum

# tag families as either gram negative or gram positive (extract taxid from description, input into taxallnomy, get phylums, list phylums in each group?)
## phylum groups got from prev step

# any gram positives in the tagged sequences at each percent confident level?


# If time:
# any other regions to tag and examine? eg low conservation across all sequences, or other 'dodgy' non-conserved regions?
## there is another region near start, shorter region

# better cleaning
# Systematic approach: for positions matching criteria, identify which sequences that are not gaps at that position and should be removed, repeat for all positions in identified conserved regions
from Bio import AlignIO
alignment = AlignIO.read("output_cleaned_v3.fa", "fasta")

new = []
for seq in alignment:
    if seq[231] != '-':
        new.append(seq.description)

print(new)
"""