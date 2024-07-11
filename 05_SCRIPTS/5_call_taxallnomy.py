import requests
from tqdm import tqdm

# read in taxonomic ids from the combined list
with open ("./01_DATA/Amidase_3/03_2_Sequence_Annotation/taxallnomy_input.txt", "r") as tax:
    tax_ids_all = tax.read().split(", ")

# API call to Taxallnomy for family ID to help cluster groups
familyID = []
print ("Now calling Taxallnomy...")

for id in tqdm(tax_ids_all):
    call = requests.get(f"http://bioinfo.icb.ufmg.br/cgi-bin/taxallnomy/taxallnomy_multi.pl?txid={id}&rank=main&format=json")

    try:
        call.status_code == 200

        
    except:
        raise Exception(f"API failed to call (status code: {call.status_code})")
        continue
    
    # pull out phylum, species, tax_ID as new unique ID
    try:
        call_json = call.json()
        new_id = (">" + f"{id}" + ":" + (call_json[f"{id}"]["phylum"]) + ":" + (call_json[f"{id}"]["family"]) + ":" + (call_json[f"{id}"]["species"]) + ",")
    except:
        new_id = (">" + f"Nil, {id}" + ",")
    
    # saves call to text file in case of issues later in workflow
    with open ("./01_DATA/Amidase_3/03_2_Sequence_Annotation/familyID_temp.txt", "a") as addfam:
        addfam.write(new_id)

# load in family ID from txt file
with open ("./01_DATA/Amidase_3/03_2_Sequence_Annotation/familyID_temp.txt", "r") as f:
    familyID = f.read().split(",>")

# replace sequenceID in combined fasta with new unique ID as above
counter = 0
print("Now appending new taxon identifier to sequences...")
with open ("./01_DATA/Amidase_3/03_1_Sequence_Searches/BLAST_Interpro_all_seq.txt", "r") as f:
    for line in f:
        if line.startswith(">"):
            line = ">" + familyID[counter] + "\n"
            counter += 1
        with open ("./01_DATA/Amidase_3/03_2_Sequence_Annotation/all_seq_taxid.txt", "a") as add:
            add.write(line) 

count = 0
with open("./01_DATA/Amidase_3/03_2_Sequence_Annotation/all_seq_taxid.txt", 'r') as file:
    for line in file:
        if line.startswith(">"):
            count += 1
print(count)