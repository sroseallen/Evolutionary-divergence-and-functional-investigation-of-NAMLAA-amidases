from Bio import SeqIO
import os

input = []
output = "all_seq.txt"

# identify all blastp sequence outputs
for file_name in os.listdir("./sequences"):
    if file_name.endswith(".txt"):
        input.append(os.path.join("./sequences", file_name))

# combine files
with open ("all_seq.txt", "w") as combo:
    for file_name in input:
        with open (file_name, "r") as f:
            combo.write(f.read())

with open ("all_seq.txt", "r") as f:
    spaceloss = f.read().replace(" ", "_")
    
with open ("all_seq_nospace.txt", "w") as f:
    f.write(spaceloss)

# remove duplicate sequences from combined file
unique_seq = {}

for record in SeqIO.parse("all_seq_nospace.txt", "fasta"):
    unique_seq[record.id] = record.seq
print(len(unique_seq))

# put the unique dictionary of sequences back in a fasta file for MSA
all = []
for idrec, seq in unique_seq.items():
    record = SeqIO.SeqRecord(seq, id = idrec, description = "")
    all.append(record)
with open ("all_seq_unique.txt", "w") as output:
    SeqIO.write (all, output, "fasta")
