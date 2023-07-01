from Bio import SeqIO
import requests
from tqdm import tqdm

full_names = []
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
    call = requests.get(f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon_suggest/{name}").json()
    try:
        taxid.append(call["sci_name_and_ids"][0]["tax_id"])
    except:
        taxid.append(f"{name} error")

# API call to Taxallnomy to cluster groups

familyID = []
print ("Now calling Taxallnomy...")
for id in tqdm(taxid):
    call = requests.get(f"http://bioinfo.icb.ufmg.br/cgi-bin/taxallnomy/taxallnomy_multi.pl?txid={id}&rank=main&format=json").json()

    # pull out phylum, species, tax_ID as new unique ID
    try:
        new_id = (">" + f"{id}" + ":" + (call[f"{id}"]["family"]) + ":" + (call[f"{id}"]["species"]) + "\n")
    except:
        new_id = (f"{id} error")
    familyID.append(new_id)

# replace sequenceID in fasta with new unique ID as above
counter = 0
print("Now appending new taxon identifier to sequences...")
with open ("all_seq_unique.txt", "r") as f:
    for line in f:
        if line.startswith(">"):
            line = familyID[counter]
            counter += 1
        open ("all_seq_taxID.txt", "a").write(line)
            
# remove outgroups from previous alignment (regions with >99.95% identity/conservation across all other positions but gaps introduced - 0.05% gapped), then re-align

# tag sequences with >50% identity in the helix 5 region

# % of identicial/consistent class/phylum in the tagged sequences? In the untagged sequences?
# any gram positives in the tagged sequences? Identify using family (eg Enterobactericae)

# any other regions to tag and examine? eg low conservation across all sequences, or other 'dodgy' non-conserved regions?