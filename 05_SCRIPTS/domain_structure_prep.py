full_names = []
full_dict = {}

# add UniProt/NCBI accessions back to sequence name (based on taxid which should be unique)
## open sequence file and filter for the fasta titles (>)
with open("./01_DATA/Amidase_3/03_1_Sequence_Searches/BLAST_Interpro_all_seq.txt", 'r') as f:
    for line in f:
        nospace = line.replace(' ', '')
        if nospace.startswith(">"):
            clean = nospace.replace(">","")
            full_names.append(clean)

## fill in NCBI taxids from the NCBI accession
from Bio import Entrez
import json
Entrez.email="xxx"

# Source for accession_to_taxid_bulk functions: https://bioinformatics.stackexchange.com/questions/4586/retrieving-ncbi-taxa-ids-from-refseq-or-genbank-assembly-accession
def accession_to_taxid_bulk(acc: list[str], db='nuccore') -> dict[str, str]:
    handle = Entrez.epost(db=db, id=','.join(acc))
    record = Entrez.read(handle)
    handle = Entrez.esummary(db=db, query_key = record['QueryKey'], WebEnv=record['WebEnv'])
    result = Entrez.read(handle)
    return {i['AccessionVersion']: int(i['TaxId']) for i in result}

## separate out the taxids and accession numbers + read into dict (accession=key, taxid=value)
import re
for_bulk_ncbi = []
for i in full_names:
    v = re.split(r'[|:]', i)[0]
    if v.startswith("WP_"):
        for_bulk_ncbi.append(v)
    else:
        k = re.split(r'[|:]', i)[3]
    full_dict[k] = v

ncbi_accessions = accession_to_taxid_bulk(for_bulk_ncbi, db="protein")
ncbi_swap = {str(value): key for key, value in ncbi_accessions.items()}

full_dict_combo = {**full_dict, **ncbi_swap}

## open cluster file + separate out the taxids
cluster_keys = []
with open ("./01_DATA/Amidase_3/05_1_Region_Conservation/[0, 0, 0, 0, 0, 0, 0, 0]","r") as cluster:
    fingerprints = json.load(cluster)
    key_list = fingerprints.keys()
    for key in key_list:
        key = key.split(":")[1]
        cluster_keys.append(key)

## make a list of accession numbers in the cluster file named after the fingerprint
accessions = []
for key in cluster_keys:
    accessions.append(full_dict_combo.get(key))
print(len(accessions))

## repeat for all cluster files, should get 62 accession lists to search with

# download full sequences for each cluster based on the accessions
## NCBI download

## InterPro download

# Pre-predicted structures
## AlphaFold DB API to get model quality estimates for the existing UniProt IDs

## Download mmCIF files for structures with a good quality estimate

# New structures (NCBI accessions)
## random selection of up to 5 sequences from the cluster file
## save full FASTA sequences in a file 
### run though OpenFold, or ColabFold in batch mode as this may be faster? 
