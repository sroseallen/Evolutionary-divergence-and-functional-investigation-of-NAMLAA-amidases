from Bio import AlignIO
from matplotlib import pyplot as plt
import os

# initialise dictionary for binary fingerprints
fingerprints = {}
id = 0
for seq in AlignIO.read(
    "./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa", "fasta"
):
    unique_key = str(id) + ":" + seq.description
    fingerprints[unique_key] = []
    id += 1


# create region fingerprint for each sequence
def region_fingerprint(
    start,
    end,
    region_length,
    average,
    alignment_path="./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa",
):
    alignment = AlignIO.read(alignment_path, "fasta")
    insertion_region = list(range(start, end))
    id = 0

    for seq in alignment:
        seq_occupancy = 0
        unique_key = str(id) + ":" + seq.description
        # quantify the number of occuped (non-gap) columns in the region
        for pos in insertion_region:
            if seq[pos - 1] != "-":
                seq_occupancy += 1
        # if occupancy more than average, consider species to have the insertion, add to fingerprint
        if seq_occupancy >= average:
            # sufficient occupancy in alignment: insertion region present
            fingerprints[unique_key].append(1)
        elif seq_occupancy < average:
            # not enough occupancy: no insertion region
            fingerprints[unique_key].append(0)
        id += 1


print("Creating fingerprinting dictionary...")
region_fingerprint(
    15,
    57,
    43,
    12.9,
    alignment_path="./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa",
)  # I-2
region_fingerprint(
    96,
    128,
    33,
    5.3,
    alignment_path="./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa",
)  # I-4
region_fingerprint(
    154,
    175,
    22,
    6.1,
    alignment_path="./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa",
)  # I-5
region_fingerprint(
    184,
    241,
    58,
    24.0,
    alignment_path="./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa",
)  # I-6
region_fingerprint(
    252,
    303,
    52,
    25.4,
    alignment_path="./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa",
)  # I-7

# this should give a dictionary of 20304 key-value pairs, where each key is the species and each value is a list of 5 0/1 values
# where 0 indicates no insertion, and 1 indicates the insertion is present

import json

with open(
    "./01_DATA/Amidase_3/05_3_Region_Conservation_24567/fingerprints.json", "w"
) as file:
    json.dump(fingerprints, file)
print("Fingerprints saved!")

# cluster species into groups based on having the same binary fingerprint
## re-reading the dictionary
import json

with open(
    "./01_DATA/Amidase_3/05_3_Region_Conservation_24567/fingerprints.json", "r"
) as file:
    fingerprints = json.load(file)

# generate filtered dictionaries for each possible cluster of insertions
import itertools

length = 5
options = list(itertools.product([0, 1], repeat=length))
options_list = [list(option) for option in options]
for option in options_list:
    cluster = {key: value for key, value in fingerprints.items() if value == option}
    if len(cluster) != 0:
        with open(
            f"./01_DATA/Amidase_3/05_3_Region_Conservation_24567/{option}", "w"
        ) as f:
            json.dump(cluster, f)


# barplot representation of each fingerprint region
length_dict = {}
for file_name in os.listdir("./01_DATA/Amidase_3/05_3_Region_Conservation_24567/"):
    if file_name.startswith("["):
        with open(
            f"./01_DATA/Amidase_3/05_3_Region_Conservation_24567/{file_name}", "r"
        ) as file:
            fingerprints = json.load(file)
            print(len(fingerprints))
            length_dict[file_name] = len(fingerprints)

length_dict = {k: v for k, v in sorted(length_dict.items(), key=lambda item: item[1])}
print(length_dict)

plt.bar(
    range(len(length_dict)), list(length_dict.values()), align="center", color="navy"
)
plt.xlabel("Fingerprint ID")
plt.ylabel("No. occurrences")
plt.title("Barplot of all fingerprint IDs")
for i in range(len(length_dict)):
    plt.text(
        i + 0.05,
        list(length_dict.values())[i] + 250,
        list(length_dict.values())[i],
        ha="center",
        va="center",
        fontsize=6,
        rotation=90,
    )
plt.xticks(range(len(length_dict)), list(length_dict.keys()), rotation=90, fontsize=6)
plt.yticks(range(0, 11001, 500))
plt.margins(x=0)
plt.grid(True, color="#EEEEEE", axis="y")
plt.show()
