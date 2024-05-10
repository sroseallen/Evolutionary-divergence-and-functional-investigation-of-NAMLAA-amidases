from Bio import AlignIO
import pandas as pd
from matplotlib import pyplot as plt

alignment = AlignIO.read("./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa", "fasta")

# set boundaries for the low_conserved region
helix_predicted_at = list(range(149, 220)) # length = 72

# write sequence annotations into dataframe
df = pd.DataFrame({
    "taxid": [],
    "family": [],
    "genus_species":[],
    "amr_fullname":[],
    "gram_status":[],
    "pathogenicity":[],
    "I6_occupancy": []
})

for seq in alignment:
    seq_occupancy = 0
    seq_split = seq.description.split(":")
    # quantify the number of occuped (non-gap) columns in the region
    for pos in helix_predicted_at:
        if seq[pos-1] != '-':
            seq_occupancy += 1
    # add this as a new annotation to the species name
    seq_split.append(seq_occupancy)
    # read annotation into the initialised dataframe above
    seq_split = pd.DataFrame([seq_split], columns=df.columns)
    df = pd.concat([df, seq_split], ignore_index=True)

df.to_csv("./01_DATA/Amidase_3/05_1_MSA_Helix_Conservation/species_data.csv")

# file inspected for unknown AMR calls and annotated using NCBI Taxonomy Browser
import pandas as pd
from matplotlib import pyplot as plt

df_2 = pd.read_csv("./01_DATA/Amidase_3/05_1_MSA_Helix_Conservation/species_data.csv")

# reformat for plotting: gram staining status 
df_grouped = df_2.groupby(["I6_occupancy","gram_status"]).size().reset_index(name="count")
df_plotting = df_grouped.pivot_table(index="I6_occupancy", columns="gram_status", values='count', fill_value=0).reset_index()

# stacked bar chart for raw counts
df_plotting.plot(x='I6_occupancy', kind='bar', stacked=True)
plt.title('Conservation across I-6 region, Raw Counts')
plt.xlabel('Number of non-gap I6 columns')
plt.ylabel('Number of species')
plt.xticks(rotation=0, ha='center')
plt.show()