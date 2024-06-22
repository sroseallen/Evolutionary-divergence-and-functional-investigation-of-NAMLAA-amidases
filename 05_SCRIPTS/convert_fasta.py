from Bio import AlignIO

input = AlignIO.parse(open("01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa", "r"), "fasta")
AlignIO.write(input, open("01_DATA/Amidase_3/07_Structure_Predictions/alignment_3_thresh1.0.a3m", "w"), "fasta-2line")

