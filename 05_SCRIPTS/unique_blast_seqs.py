from Bio import SeqIO
import os

input = []
output = "./01_DATA/Amidase_3/03_1_Sequence_Searches/BLAST_all_seq.txt"

# identify all blastp sequence outputs
for file_name in os.listdir("./01_DATA/Amidase_3/03_1_Sequence_Searches/"):
    if file_name.endswith("seqdump.txt"):
        input.append(os.path.join("./01_DATA/Amidase_3/03_1_Sequence_Searches/", file_name))

# combine files
with open ("./01_DATA/Amidase_3/03_1_Sequence_Searches/BLAST_all_seq.txt", "w") as combo:
    for file_name in input:
        with open (file_name, "r") as f:
            combo.write(f.read())

with open ("./01_DATA/Amidase_3/03_1_Sequence_Searches/BLAST_all_seq.txt", "r") as f:
    spaceloss = f.read().replace(" ", "_")
    
with open ("./01_DATA/Amidase_3/03_1_Sequence_Searches/BLAST_all_seq_nospace.txt", "w") as f:
    f.write(spaceloss)

# remove duplicate sequences from combined file
unique_seq = {}

for record in SeqIO.parse("./01_DATA/Amidase_3/03_1_Sequence_Searches/BLAST_all_seq_nospace.txt", "fasta"):
    unique_seq[record.id] = record.seq
print(len(unique_seq))

# put the unique dictionary of sequences back in a fasta file for MSA
all = []
for idrec, seq in unique_seq.items():
    record = SeqIO.SeqRecord(seq, id = idrec, description = "")
    all.append(record)
with open ("./01_DATA/Amidase_3/03_1_Sequence_Searches/BLAST_all_seq_unique.txt", "w") as output:
    SeqIO.write (all, output, "fasta")
