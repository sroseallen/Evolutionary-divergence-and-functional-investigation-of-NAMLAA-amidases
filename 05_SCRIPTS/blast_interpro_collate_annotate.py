from Bio import SeqIO
import requests
from tqdm import tqdm

InterPro_taxid = "./01_DATA/Amidase_3/03_1_Sequence_Searches/InterPro_IPR002508_API_seqs.txt"
BLAST_taxid = "./01_DATA/Amidase_3/03_2_Sequence_Annotation/taxid_temp.txt"
BLAST_seqs = "./01_DATA/Amidase_3/03_1_Sequence_Searches/BLAST_all_seq_unique.txt"
all_seqs = "./01_DATA/Amidase_3/03_1_Sequence_Searches/BLAST_Interpro_all_seq.txt"

# read in BLAST taxids
with open (BLAST_taxid, "r") as blast:
    blast_ids = blast.read().split(", ")
    blast_id_set = set(blast_ids)

tax_ids = []
uniquetaxis = 0

# identify new unique sequences from InterPro download
records = list(SeqIO.parse(InterPro_taxid, "fasta"))

with open (all_seqs, "a") as allseqs:
    for record in records:
        taxid = record.description.split("|")[3]
        if taxid not in blast_id_set and taxid not in tax_ids:
            tax_ids.append(taxid)
            uniquetaxis += 1
            SeqIO.write(record, allseqs, "fasta")
print(uniquetaxis)
print(len(tax_ids))

# add all BLAST RefSeq sequences and ID lines to the all_seqs list
with open (BLAST_seqs, "r") as blast, open (all_seqs, "a") as allseqs:
    for line in blast:
        allseqs.write(line)

# create list of all taxids for this combined list
tax_ids.extend(blast_ids)
print(len(tax_ids))

# API call to Taxallnomy for family ID to help cluster groups
familyID = []
print ("Now calling Taxallnomy...")
for id in tqdm(tax_ids):
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