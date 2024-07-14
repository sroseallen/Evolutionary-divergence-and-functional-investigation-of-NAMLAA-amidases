full_names = []
full_dict = {}

# add UniProt/NCBI accessions back to sequence name (based on taxid which should be unique)
## open sequence file and filter for the fasta titles (>)
with open(
    "./01_DATA/Amidase_3/03_1_Sequence_Searches/BLAST_Interpro_all_seq.txt", "r"
) as f:
    for line in f:
        nospace = line.replace(" ", "")
        if nospace.startswith(">"):
            clean = nospace.replace(">", "")
            full_names.append(clean)

## fill in NCBI taxids from the NCBI accession
from Bio import Entrez
import json

Entrez.email = "XXX"


# Source for accession_to_taxid_bulk functions: https://bioinformatics.stackexchange.com/questions/4586/retrieving-ncbi-taxa-ids-from-refseq-or-genbank-assembly-accession
def accession_to_taxid_bulk(acc: list[str], db="nuccore") -> dict[str, str]:
    handle = Entrez.epost(db=db, id=",".join(acc))
    record = Entrez.read(handle)
    handle = Entrez.esummary(
        db=db, query_key=record["QueryKey"], WebEnv=record["WebEnv"]
    )
    result = Entrez.read(handle)
    return {i["AccessionVersion"]: int(i["TaxId"]) for i in result}


## separate out the taxids and accession numbers + read into dict (accession=key, taxid=value)
import re

for_bulk_ncbi = []
for i in full_names:
    v = re.split(r"[|:]", i)[0]
    if v.startswith("WP_"):
        for_bulk_ncbi.append(v)
    else:
        k = re.split(r"[|:]", i)[3]
    full_dict[k] = v

# read in ncbi accessions that are missing
ncbi_accessions = accession_to_taxid_bulk(for_bulk_ncbi, db="protein")
ncbi_swap = {str(value): key for key, value in ncbi_accessions.items()}

full_dict_combo = {**full_dict, **ncbi_swap}

# save accession:taxid dictionary
with open(
    "./01_DATA/Amidase_3/07_Structure_Predictions/accession_taxid_dictionary.json", "w"
) as f:
    json.dump(full_dict_combo, f)


## open cluster file + separate out the taxids
def gen_accession(filepath: str) -> list:
    cluster_keys = []
    with open(filepath, "r") as cluster:
        fingerprints = json.load(cluster)
        key_list = fingerprints.keys()
        for key in key_list:
            key = key.split(":")[1]
            cluster_keys.append(key)

    ## make a list of accession numbers in the cluster file named after the fingerprint
    accessions = []
    for key in cluster_keys:
        accessions.append(full_dict_combo.get(key))
    return accessions


## repeat for all cluster files, should get 62 accession lists to search with
import os

accession_dict = {}

for file_name in os.listdir("./01_DATA/Amidase_3/05_1_Region_Conservation/"):
    if file_name.startswith("["):
        cluster = gen_accession(
            f"./01_DATA/Amidase_3/05_1_Region_Conservation/{file_name}"
        )
        print(f"{file_name}:", len(cluster))
        accession_dict[f"{file_name}"] = cluster


# Pre-predicted structures (UniProt only)
import requests
from tqdm import tqdm

## Download mmCIF files for structures with a UniProt ID
list_fingerprints = []
for file_name in os.listdir("./01_DATA/Amidase_3/05_1_Region_Conservation/"):
    if file_name.startswith("["):
        list_fingerprints.append(file_name)


def download_strucs(accs: list, fingerprint_concat: str):
    base_url = "https://alphafold.ebi.ac.uk/api/prediction/"

    print("Now calling AlphaFold Database...")
    for name in tqdm(accs):
        try:
            if not name.startswith("WP_"):
                call = requests.get(f"{base_url}{name}")
                data = call.json()
                # download mmCIF file from AlphaFold Database
                name = f"./01_DATA/Amidase_3/07_Structure_Predictions/{fingerprint_concat}_{name}.cif"
                response = requests.get(data[0]["cifUrl"])
                if response.status_code == 200:
                    with open(name, "wb") as f:
                        f.write(response.content)
            else:
                continue

        except AssertionError as e:
            print("Assertion Error has happened - issue with accession ID")
            continue

        except AttributeError as a:
            print("Assertion Error has happened - issue with accession ID")
            continue

        except:
            print(
                f"API failed, {name} not in database (status code: {call.status_code})"
            )
            continue


for i in list_fingerprints:
    accs = accession_dict.get(i)
    fingerprint_concat = "".join(re.findall(r"\d+", i))
    download_strucs(accs=accs, fingerprint_concat=fingerprint_concat)

# New structures (process)
## download full sequences for 5 random sequences from top 7n clusters:
import random


def select_seqs(key: str) -> list:
    valid_accessions = [
        value for value in list(accession_dict.get(key)) if value != None
    ]
    selection = random.sample(valid_accessions, 5)
    return selection


rep_seqs = []

# 00001111
rep_seqs.append(select_seqs("[0, 0, 0, 0, 1, 1, 1, 1]"))
# 00001110
rep_seqs.append(select_seqs("[0, 0, 0, 0, 1, 1, 1, 0]"))
# 00001000
rep_seqs.append(select_seqs("[0, 0, 0, 0, 1, 0, 0, 0]"))
# 01010000
rep_seqs.append(select_seqs("[0, 1, 0, 1, 0, 0, 0, 0]"))
# 00011000
rep_seqs.append(select_seqs("[0, 0, 0, 1, 1, 0, 0, 0]"))
# 00000000
rep_seqs.append(select_seqs("[0, 0, 0, 0, 0, 0, 0, 0]"))
# 00100000
rep_seqs.append(select_seqs("[0, 0, 1, 0, 0, 0, 0, 0]"))
# 00101110
rep_seqs.append(select_seqs("[0, 0, 1, 0, 1, 1, 1, 0]"))
# 00101000
rep_seqs.append(select_seqs("[0, 0, 1, 0, 1, 0, 0, 0]"))

# 00011100 - I4 and I6
rep_seqs.append(select_seqs("[0, 0, 0, 1, 1, 1, 0, 0]"))
# 10001111 - I1 and the I5/I6/I7 combo
rep_seqs.append(select_seqs("[1, 0, 0, 0, 0, 1, 1, 0]"))
# 00100111 - I6/I7 without I5
rep_seqs.append(select_seqs("[0, 0, 1, 0, 0, 1, 1, 1]"))

## save accessions in a file (to search and get full sequences for predictions)
with open(
    "./01_DATA/Amidase_3/07_Structure_Predictions/accessions_to_search.txt", "a"
) as f:
    f.write(",".join(str(i) for i in rep_seqs))
