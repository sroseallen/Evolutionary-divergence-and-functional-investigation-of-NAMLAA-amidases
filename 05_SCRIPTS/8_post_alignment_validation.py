from Bio import SeqIO
import pandas as pd
import seaborn as sn 
import matplotlib.pyplot as plt 

seqs = list(SeqIO.parse("./01_DATA/Amidase_3/04_Multiple_Alignment/Validation/candidate_validation_set.txt", "fasta"))
strucs = list(SeqIO.parse("./01_DATA/Amidase_3/04_Multiple_Alignment/Validation/candidate_USalign_validation_set.txt", "fasta-2line"))

structures = ["4BIN_","3NE8_","4RN7_","7RAG_","7TJ4_","5EMI_","4M6G_","4LQ6_","1XOV_"]

# compare sequences and structures in a pairwise fashion
def validate_msa(A, B):
    for record in seqs:
        if A in record.id: msa1_temp = record.seq
        if B in record.id: msa2_temp = record.seq

    # msa: remove gaps aligned with gaps for seq1 and seq2
    msa1 = []
    msa2 = []
    # count gaps that will remain after removing gap-gap pairs
    msa_gaps = 0
    # index = base index (gaps have an index number)
    index = 0
    # index_counter = indices with information (gaps do not add to the index_counter)
    index_counter = 0
    aligned_index_seq = []

    for char in msa1_temp:
        if char != "-": 
            msa1.append(char)
            index_counter += 1
            # if char is not a gap, check if its aligned to a gap. If it isn't, then the seq pair is 'aligned'.
            if msa2_temp[index] != "-":
                # add index to index list,then move to next char in seq1
                aligned_index_seq.append(index_counter)
            # if char is aligned to a gap, them do not add to the index list.
            elif msa2_temp[index] == "-":
                next
        # if char itself is a gap, do not add to aligned index list or to the index counter, but still add to the msa1 list if its aligned to a non-gap.
        elif char == "-":
            if msa2_temp[index] != "-":
                msa1.append(char)
                msa_gaps += 1
        index += 1

    index_counter=0
    for char in msa2_temp:
        if char != "-":
            msa2.append(char)
        elif char == "-":
            if msa1_temp[index_counter] != "-":
                msa2.append(char)
                msa_gaps += 1
        index_counter += 1

    # total aligned pairs in msa
    msa_aligned_pairs = len(msa1) - msa_gaps

    # length of alignment between 4BIN and 3NE8 in structural alignment
    for record in strucs:
        if record.id == str(A+B+"USalign"): struc1 = record.seq
        if record.id == str(B+A+"USalign"): struc2 = record.seq

    # count number of gaps "-" in seq1 and seq2 in structural alignment, add together.
    struc_gaps=0
    # index = base index (gaps have an index number)
    index_struc = 0
    # index_counter = indixes with information (gaps do not add to the index_counter)
    index_counter_struc = 0
    aligned_index_struc = []

    for char in struc1:
        # because 4BIN structure is incomplete:
        if index_counter_struc == 120 and A == "4BIN_":
            index_counter_struc+=12
        if index_counter_struc == 121 and A == "3NE8_":
            index_counter_struc+=6
        if char != "-": 
            index_counter_struc += 1
            # if char is not a gap, check if its aligned to a gap. If it isn't, then the seq pair is 'aligned'.
            if struc2[index_struc] != "-":
                # add index to index list,then move to next char in seq1
                aligned_index_struc.append(index_counter_struc)
            # if char is aligned to a gap, them do not add to the index list.
            elif struc2[index_struc] == "-":
                next
        # if char itself is a gap, do not add to aligned index list or to the index counter
        elif char == "-":
            struc_gaps+=1
        index_struc += 1

    for char in struc2:
        if char == "-":
            struc_gaps+=1

    # total aligned pairs in structural alignment
    struc_aligned_pairs = len(struc1) - struc_gaps

    # number of pairs aligned in msa and in struc alignment
    aligned_seq_struct = 0
    for indice in aligned_index_seq:
        if indice in aligned_index_struc:
            aligned_seq_struct+=1

    # calculate TPR and PPV
    TPR = aligned_seq_struct/struc_aligned_pairs
    PPV = aligned_seq_struct/msa_aligned_pairs
    return[TPR, PPV]

# initialise dataframes for storing TPR and PPV values for each pair
tpr = pd.DataFrame(
    index=structures,
    columns=structures
)
ppv = pd.DataFrame(
    index=structures,
    columns=structures
)

# generate values for each pair
for i,s in enumerate(structures):
    A = structures[i]
    for j,t in enumerate(structures):
        B = structures[j]
        tpr.loc[A,B]  = float(validate_msa(A,B)[0])
        ppv.loc[A,B]  = validate_msa(A,B)[1]

for col in structures:
    tpr[col]=tpr[col].astype(float)
    ppv[col]=ppv[col].astype(float)

# remove redundant top triangle of results in the heatmap plot (tpr plot)
import numpy as np
mask = np.zeros_like(tpr, dtype=bool)
mask[np.triu_indices_from(mask)] = True
mask[np.diag_indices_from(mask)] = False

# save tpr heatmaps
print(tpr.values.mean())
sn.heatmap(data=tpr,square=True,annot=True,fmt='.2f',mask=mask,cmap="crest",vmin=0.8,vmax=1.0)
plt.tight_layout()
plt.savefig("./01_DATA/Amidase_3/04_Multiple_Alignment/Validation/tpr.png")
plt.close()

# repeat the above for the ppv plot
mask = np.zeros_like(ppv, dtype=bool)
mask[np.triu_indices_from(mask)] = True
mask[np.diag_indices_from(mask)] = False

print(ppv.values.mean())
sn.heatmap(data=ppv,square=True,annot=True,fmt='.2f',mask=mask,cmap="crest",vmin=0.8,vmax=1.0)
plt.tight_layout()
plt.savefig("./01_DATA/Amidase_3/04_Multiple_Alignment/Validation/ppv.png")