from Bio import AlignIO
import pandas as pd
from matplotlib import pyplot as plt

def plot_insertion(start, end, save_path, alignment_path = "./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa", write=False):
    alignment = AlignIO.read(alignment_path, "fasta")
    insertion_region = list(range(start, end))

    # write sequence annotations into dataframe
    df = pd.DataFrame({
        "taxid": [],
        "phylum": [],
        "family": [],
        "genus_species":[],
        "gram_status":[],
        "occupancy": []
    })

    for seq in alignment:
        seq_occupancy = 0
        seq_split = seq.description.split(":")
        # quantify the number of occuped (non-gap) columns in the region
        for pos in insertion_region:
            if seq[pos-1] != '-':
                seq_occupancy += 1
        # add this as a new annotation to the species name
        seq_split.append(int(seq_occupancy))
        # read annotation into the initialised dataframe above
        seq_split = pd.DataFrame([seq_split], columns=df.columns)
        df = pd.concat([df, seq_split], ignore_index=True)
    
    if write==True:
        df.to_csv(f"./01_DATA/Amidase_3/05_1_Region_Conservation/{start}-{end}.csv", index=False)
            
    # reformat for plotting: gram staining status 
    df_grouped = df.groupby(["occupancy","gram_status"]).size().reset_index(name="count")
    df_grouped["occupancy"] = df_grouped["occupancy"].astype("int64")
    for num in range((end-start)+2):
        df_grouped = pd.concat([df_grouped, pd.DataFrame([{"occupancy":num, "gram_status":"NA", "count":0}])])
    df_plotting = df_grouped.pivot_table(index="occupancy", columns="gram_status", values='count', fill_value=0).reset_index()

    # stacked bar chart for raw counts
    df_plotting.plot(x='occupancy', kind='bar', stacked=True)
    plt.title('Conservation across region, Raw Counts')
    plt.xlabel('Number of non-gap columns')
    plt.ylabel('Number of species')
    plt.xticks(rotation=90, ha='center')
    plt.tight_layout()
    plt.savefig(save_path)

print ("Now plotting region 1...")
plot_insertion(1, 6, save_path="./01_DATA/Amidase_3/05_1_Region_Conservation/I1_region.png") #I-1
print ("Now plotting region 2...")
plot_insertion(15, 57, save_path="./01_DATA/Amidase_3/05_1_Region_Conservation/I2_region.png") #I-2
print ("Now plotting region 3...")
plot_insertion(81, 85, save_path="./01_DATA/Amidase_3/05_1_Region_Conservation/I3_region.png") #I-3
print ("Now plotting region 4...")
plot_insertion(96, 128, save_path="./01_DATA/Amidase_3/05_1_Region_Conservation/I4_region.png") #I-4
print ("Now plotting region 5...")
plot_insertion(154, 175, save_path="./01_DATA/Amidase_3/05_1_Region_Conservation/I5_region.png") #I-5
print ("Now plotting region 6...")
plot_insertion(184, 241, save_path="./01_DATA/Amidase_3/05_1_Region_Conservation/I6_region.png",write=True) #I-6
print ("Now plotting region 7...")
plot_insertion(252, 303, save_path="./01_DATA/Amidase_3/05_1_Region_Conservation/I7_region.png",write=True) #I-7
print ("Now plotting region 8...")
plot_insertion(355, 399, save_path="./01_DATA/Amidase_3/05_1_Region_Conservation/I8_region.png") #I-8
print ("Plots saved!")


from Bio import AlignIO
import pandas as pd
from matplotlib import pyplot as plt

# initialise dictionary for binary fingerprints
fingerprints = {}
id=0
for seq in AlignIO.read("./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa", "fasta"):
    unique_key = str(id) + ":" + seq.description
    fingerprints[unique_key] = []
    id +=1

# create region fingerprint for each sequence
def region_fingerprint(start, end, region_length, average, alignment_path = "./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa"):
    alignment = AlignIO.read(alignment_path, "fasta")
    insertion_region = list(range(start, end))
    id=0

    for seq in alignment:
        seq_occupancy = 0
        unique_key = str(id) + ":" + seq.description
        # quantify the number of occuped (non-gap) columns in the region
        for pos in insertion_region:
            if seq[pos-1] != '-':
                seq_occupancy += 1
        # if occupancy more than average, consider species to have the insertion, add to fingerprint
        if seq_occupancy >= average:
            #sufficient occupancy in alignment: insertion region present
            fingerprints[unique_key].append(1)
        elif seq_occupancy < average:
            #not enough occupancy: no insertion region
            fingerprints[unique_key].append(0)
        id+=1

print ("Creating fingerprinting dictionary...")
region_fingerprint(1, 6, 6, 1.17, alignment_path = "./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa") #I-1
region_fingerprint(15, 57, 43, 12.9, alignment_path = "./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa") #I-2
region_fingerprint(81, 85, 5, 1.3, alignment_path = "./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa") #I-3
region_fingerprint(96, 128, 33, 5.3, alignment_path = "./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa") #I-4
region_fingerprint(154, 175, 22, 6.1, alignment_path = "./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa") #I-5
region_fingerprint(184, 241, 58, 24.0, alignment_path = "./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa") #I-6
region_fingerprint(252, 303, 52, 25.4, alignment_path = "./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa") #I-7
region_fingerprint(355, 399, 45, 6.1, alignment_path = "./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa") #I-8

# this should give a dictionary of 20304 key-value pairs, where each key is the species and each value is a list of 8 0/1 values 
# where 0 indicates no insertion, and 1 indicates the insertion is present

import json
with open ("./01_DATA/Amidase_3/05_1_Region_Conservation/fingerprints.json", "w") as file:
    json.dump(fingerprints, file)
print ("Fingerprints saved!")

# cluster species into groups based on having the same binary fingerprint
## re-reading the dictionary
import json
with open("./01_DATA/Amidase_3/05_1_Region_Conservation/fingerprints.json", "r") as file:
    fingerprints = json.load(file)

# generate filtered dictionaries for each possible cluster of insertions
import itertools
length = 8
options = list(itertools.product([0, 1], repeat=length))
options_list = [list(option) for option in options]
for option in options_list:
    cluster = {key: value for key, value in fingerprints.items() if value == option}
    if len(cluster) != 0:
        with open (f"./01_DATA/Amidase_3/05_1_Region_Conservation/{option}", "w") as f:
            json.dump(cluster, f)

# randomly select n=5 sequences from each cluster
import random
import json
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

# open each fingerprint file
def fingerprint_selection (alignment_path: str) -> list:
    random_selection = {}
    with open(f"{alignment_path}", "r") as file:
        fingerprints = json.load(file)
        # obtain fingerprint as a string
        value = str(next(iter(fingerprints.values())))
        binary_print = ''.join(value.split(", "))
        binary_print = binary_print.replace("[", "").replace("]", "")

        # random selection (10 keys)
        if len(fingerprints) < 10:
            selection = fingerprints.keys()
        elif len(fingerprints) >= 10:
            selection = random.sample(list(fingerprints.keys()), 10)
        
        # change key name to note taxid and fingerprint
        for item in selection:
            id = item.split(":")[1]
            random_selection[id] = binary_print
    return random_selection

# loop through each fingerprint file
taxid_lookup = {}
import os
for file_name in os.listdir("./01_DATA/Amidase_3/05_1_Region_Conservation/"):
    if file_name.startswith("["):
        random_selection = fingerprint_selection(f"./01_DATA/Amidase_3/05_1_Region_Conservation/{file_name}")
        taxid_lookup.update(random_selection)

# save list of taxids in case needed for references
with open ("./01_DATA/Amidase_3/06_Phylogeny/representative_ids.txt", "w") as file:
    file.write(','.join(taxid_lookup))

# look up the sequence using the taxid in the MSA file alignment_3_thresh1.0
alignment=AlignIO.read("./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa","fasta")

for id in taxid_lookup.keys():
    for record in alignment:
        if record.id.split(":")[0] == id:
            record.id = f"{taxid_lookup[id]}:{record.id}"

filtered_alignment = [record for record in alignment if record.id.split(":")[1] in list(taxid_lookup.keys())]

# remove duplicate taxid/fingerprint combinations
unique_seqs = []
filtered_alignment_nodup = {}

for record in filtered_alignment:
    # unique combinations with fingerprint [0] and taxid [1]
    name = f"{record.id.split(":")[0]}:{record.id.split(":")[1]}"
    if name not in unique_seqs:
        filtered_alignment_nodup[record.id] = record
        unique_seqs.append(name)

filtered_alignment_nodup = MultipleSeqAlignment(filtered_alignment_nodup.values())
AlignIO.write(filtered_alignment_nodup, "./01_DATA/Amidase_3/06_Phylogeny/representative_alignment.fasta", "fasta")

# prepare the sequences for phylogeny (phylip format)
filtered_alignment_nodup = AlignIO.read("./01_DATA/Amidase_3/06_Phylogeny/representative_alignment.fasta", "fasta")
unique_names = set()

for record in filtered_alignment_nodup:
    # cut the name down to size for the phylip format (10 characters)
    truncated_name = f"{record.id.split(":")[0]}1"

    # make the name unique
    while truncated_name in unique_names:
        current_number = int(truncated_name[8:])
        new_number = current_number + 1
        truncated_name = f"{truncated_name[:8]}{new_number}"

    unique_names.add(truncated_name)
    record.id=truncated_name
    print(truncated_name)

# Write the new alignment
AlignIO.write(filtered_alignment_nodup, "./01_DATA/Amidase_3/06_Phylogeny/representative_phylip.phy", "phylip")
