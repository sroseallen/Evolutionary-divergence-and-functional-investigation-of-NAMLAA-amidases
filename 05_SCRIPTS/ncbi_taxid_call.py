from Bio import SeqIO
import requests
from tqdm import tqdm
# Documentation (NCBI REST API v2: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/rest-api/)

full_names = []
genus_species = []

# read in all sequence names from fasta file
with open ("BLAST_all_seq_unique.txt", "r") as f:
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

    try:
        call = requests.get(f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{name}")
    except:
        raise Exception(f"API failed to call (status code: {call.status_code})")

    try:
        call_json = call.json()
        id = str(call_json["taxonomy_nodes"][0]["taxonomy"]["tax_id"])
        taxid.append(str(id))
    except:
        id = f"{name}_error"
        taxid.append(id)

    # iteratively saves tax_ids to text file in case of loop break/internet issue while running
    with open ("taxid_temp.txt", "a") as addtax:
        addtax.write(id + ", ")