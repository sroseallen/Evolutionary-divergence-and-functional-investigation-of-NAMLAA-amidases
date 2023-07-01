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
for name in tqdm(genus_species):
    call = requests.get(f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon_suggest/{name}").json()
    try:
        taxid.append(call["sci_name_and_ids"][0]["tax_id"])
    except:
        taxid.append(f"{name} error")
print(taxid[1:10])

