from Bio import SeqIO
from Bio import AlignIO
import requests
from tqdm import tqdm
import pandas as pd

""" full_names = []
genus_species = []

# read in all sequence names from fasta file
with open ("all_seq_unique.txt", "r") as f:
    for sequence in SeqIO.parse(f, "fasta"):
        full_names.append(sequence.id)

# pull just genus and species name from fasta identifier
for name in full_names:
    index = (name.index("[")) + 1
    split = name[index:-1].replace("_", " ")
    genus_species.append(split)

# Use genus, species to API call NCBI Datasets v2 REST API and retrieve tax_id for that species.
taxid = []
print ("Now calling NCBI server...")
for name in tqdm(genus_species):

    call = requests.get(f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon_suggest/{name}")

    try:
        call.status_code == 200  
    except:
        raise Exception(f"API failed to call (status code: {call.status_code})")
        continue

    try:
        call_json = call.json()
        taxid.append(call_json["sci_name_and_ids"][0]["tax_id"])
    except:
        taxid.append(f"{name} error")

    # saves call to text file in case of issues later in workflow
    with open ("taxid_temp.txt", "a") as addtax:
        addtax.write(taxid)

# load in taxid from txt file
with open ("taxid_temp.txt", "r") as f:
    taxid = f.read().split(", ")

# API call to Taxallnomy to cluster groups
# familyID = []
print ("Now calling Taxallnomy...")
for id in tqdm(taxid):
    call = requests.get(f"http://bioinfo.icb.ufmg.br/cgi-bin/taxallnomy/taxallnomy_multi.pl?txid={id}&rank=main&format=json")

    try:
        call.status_code == 200
    except:
        raise Exception(f"API failed to call (status code: {call.status_code})")
        continue
    
    # pull out phylum, species, tax_ID as new unique ID
    try:
        call_json = call.json()
        new_id = (">" + f"{id}" + ":" + (call_json[f"{id}"]["family"]) + ":" + (call_json[f"{id}"]["species"]) + ",")
    except:
        new_id = (">" + f"Nil, {id}" + ",")
    
    # saves call to text file in case of issues later in workflow
    with open ("familyID_temp.txt", "a") as addfam:
        addfam.write(new_id)

# load in family ID from txt file
with open ("familyID_temp.txt", "r") as f:
    familyID = f.read().split(",>")

# replace sequenceID in fasta with new unique ID as above
counter = 0
print("Now appending new taxon identifier to sequences...")
with open ("all_seq_unique.txt", "r") as f:
    for line in f:
        if line.startswith(">"):
            line = ">" + familyID[counter] + "\n"
            counter += 1
        with open ("all_seq_taxID.txt", "a") as add:
            add.write(line) """

# 383 errors (genus/species tag not sufficient for taxid api call) - removed (N = 34361)

# Improve alignment: remove outliers from previous alignment (regions with <0.01% occupancy in regions of otherwise >99.9% occupancy and high consensus/identity), then re-align
## Outliers removed due to above rule: 
# WP_183078710.1 Sphingobacterium puteale (x1)
# WP_224610223.1 Deinococcus multiflagellatus (x1)
# WP_093795350.1 Sporomusa acidovorans (x2 removed)
# WP_043729524.1 Nocardia_asiatica partial seq. (x1)
# WP_225220617.1:3-132_N-acetylmuramoyl-L-alanine_amidase_[Oerskovia_merdavium] (x2 removed)
# WP_211270272.1:1-94_N-acetylmuramoyl-L-alanine_amidase,_partial_[Rhodococcus_phenolicus] (x2 removed)
# WP_141015564.1:64-157_N-acetylmuramoyl-L-alanine_amidase_[Nocardioides_sambongensis] (x2 removed)
# WP_214296915.1:203-396_N-acetylmuramoyl-L-alanine_amidase_[Geobacter_chapellei] (x3 removed)
# WP_054689895.1:3-187_N-acetylmuramoyl-L-alanine_amidase_[Desulfosarcina_cetonica] (x1)
# WP_167457244.1:317-476_N-acetylmuramoyl-L-alanine_amidase_[Sphingobacterium_detergens] (x1)
# WP_166232615.1:36-178_N-acetylmuramoyl-L-alanine_amidase_[Propioniciclava_coleopterorum] (x1)

# also removed: Archea, Eukaryota, Virus (non-bacterial kingdom sequences)
## WP_048045911.1:36-222_N-acetylmuramoyl-L-alanine_amidase_[Methanosarcina_mazei] (archea that produces methane, removed)
## Pseudococcidae x4 (mealybug, insect - cells contain endosymbiotic bacteria, removed as unclear bacterial origin)
## Corticovirus PM2 x15 (bacteriophages, not removed as relevant to project)
## Keylargovirus JL001 (bacteriophage family, not removed)
## Thermococcaceae x29 (archea, extreme thermophiles, removed)

# New Sequence total N = 34306

# re-align with kalign and see if alignment any better in conserved regions. 
# Second pass:
# removed all partial sequences
# 264641:Desulfitibacteraceae:Desulfitibacter LEGKKIYIDPGHGGTYNEGIARCGTYQVGATGQYTGQMEKNTALDIAFWLREYLKGALATVYMSRINDSAVCLWERTQEANSLNADIFVSLHHNGADDPNVSGASTHWYKSVDKPLAEVVLNSLLANTEFSPWGTGLRYNNYHVLRETNMPAVLIENGFFSNEQDDRYVYGNGKDFGNSPINYNRRDIAYAIYWGIWNYF
# 2764329:Peptostreptococcaceae:Paeniclostridium hominis
# 215743:Roseobacteraceae:Roseovarius mucosus
# and others (see output_cleaned_v3, edited in Jalview)

# New sequence total N = 34008

# tag sequences with identity in the helix 5 region
# Based on predicted helix region from 3NE8 the sequence pattern is DAIAKSLAESENKVDLLDG.....DILLDLTRRET. treat these as helical positions, anything that aligns in >90% is high confirmation, >70% as medium confirmation, >50% as low confirmation
# In current alignment: Helix 1: 754, 757, 760, 767, 769, 771, 775, 777, 780, 789, 791, 792, 799, 803, 806, 810, 815, 821, 823. Helix 2: 897, 900, 901, 902, 903, 905, 907, 908, 911, 913, 917

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

#THURSDAY
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
""" from Bio import AlignIO
alignment = AlignIO.read("output_cleaned_v3.fa", "fasta")

new = []
for seq in alignment:
    if seq[231] != '-':
        new.append(seq.description)

print(new)

 """