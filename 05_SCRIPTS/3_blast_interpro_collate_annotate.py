from Bio import SeqIO

InterPro_taxid = (
    "./01_DATA/Amidase_3/03_1_Sequence_Searches/InterPro_IPR002508_API_seqs.txt"
)
BLAST_taxid = "./01_DATA/Amidase_3/03_2_Sequence_Annotation/taxid_temp.txt"
BLAST_seqs = "./01_DATA/Amidase_3/03_1_Sequence_Searches/BLAST_all_seq_unique.txt"
all_seqs = "./01_DATA/Amidase_3/03_1_Sequence_Searches/BLAST_Interpro_all_seq.txt"

# read in BLAST taxids
with open(BLAST_taxid, "r") as blast:
    blast_ids = blast.read().split(", ")
    blast_id_set = set(blast_ids)

tax_ids = []
uniquetaxis = 0

# identify new unique sequences from InterPro download
records = list(SeqIO.parse(InterPro_taxid, "fasta"))

with open(all_seqs, "a") as allseqs:
    for record in records:
        taxid = record.description.split("|")[3]
        if taxid not in blast_id_set and taxid not in tax_ids:
            tax_ids.append(taxid)
            uniquetaxis += 1
            SeqIO.write(record, allseqs, "fasta")
print(uniquetaxis)
print(len(tax_ids))

# add all BLAST RefSeq sequences and ID lines to the all_seqs list
with open(BLAST_seqs, "r") as blast, open(all_seqs, "a") as allseqs:
    for line in blast:
        allseqs.write(line)

# create list of all taxids for this combined list
with open(
    "./01_DATA/Amidase_3/03_2_Sequence_Annotation/taxid_temp.txt", "r"
) as blast, open(
    "./01_DATA/Amidase_3/03_2_Sequence_Annotation/taxallnomy_input.txt", "w"
) as taxallnomy_input:
    blast_alltax = blast.read().split(", ")
    tax_ids_all = tax_ids + blast_alltax
    print(len(tax_ids_all))
    taxallnomy_input.write(",".join(tax_ids_all))
