from Bio import SeqIO
#import pandas as pd

# annotated seqeunces
all_seqs = "./01_DATA/Amidase_3/03_2_Sequence_Annotation/all_seq_taxid.txt"
records = list(SeqIO.parse(all_seqs, "fasta"))
keep = []
remove = []

# replace call errors with manually identified species/family/taxid information
for record in records:
    if record.description == "Nil, Pseudoclostridium thermosuccinogenes_error": record.description = ">84032:Clostridiaceae:Clostridium thermosuccinogenes"
    if record.description == "Nil, Halomonas taeanensis_error": record.description = ">284577:Halomonadaceae:Onishia taeanensis"
    if record.description == "Nil, Halomonas montanilacus_error": record.description = ">2282305:Halomonadaceae:Billgrantia montanilacus"
    if record.description == "Nil, Halomonas nanhaiensis_error": record.description = ">1258546:Halomonadaceae:Vreelandella nanhaiensis"
    if record.description == "Nil, Halomonas nigrificans_error": record.description = ">2042704:Halomonadaceae:Vreelandella nigrificans"
    if record.description == "Nil, Halomonas zhanjiangensis_error": record.description = ">1121960:Halomonadaceae:Vreelandella zhanjiangensis"
    if record.description == "Nil, Caballeronia concitans_error": record.description = ">1777133:Burkholderiaceae:Caballeronia concitans"
    if record.description == "Nil, Halomonas utahensis_error": record.description = ">86177:Halomonadaceae:Vreelandella utahensis"
    if record.description == "Nil, Halomonas aquamarina_error": record.description = ">77097:Halomonadaceae:Vreelandella aquamarina"
    if record.description == "Nil, Halomonas malpeensis_error": record.description = ">1172368:Halomonadaceae:Vreelandella malpeensis"
    if record.description == "Nil, Halomonas gudaonensis_error": record.description = ">376427:Halomonadaceae:Billgrantia gudaonensis"
    if record.description == "Nil, Halomonas hamiltonii_error": record.description = ">502829:Halomonadaceae:Vreelandella hamiltonii"
    if record.description == "Nil, Halomonas maris_error": record.description = ">2729617:Halomonadaceae:Vreelandella maris"
    if record.description == "Nil, Halomonas andesensis_error": record.description = ">447567:Halomonadaceae:Vreelandella andesensis"
    if record.description == "Nil, Halomonas kenyensis_error": record.description = ">321266:Halomonadaceae:Billgrantia kenyensis"
    if record.description == "Nil, Halomonas populi_error": record.description = ">2498858:Halomonadaceae:Vreelandella populi"
    if record.description == "Nil, Halomonas sedimenti_error": record.description = ">2729618:Halomonadaceae:Vreelandella sedimenti"
    if record.description == "Nil, Halomonas sulfidaeris_error": record.description = ">115553:Halomonadaceae:Vreelandella sulfidaeris"
    if record.description == "Nil, Halomonas alkaliphila_error": record.description = ">272774:Halomonadaceae:Vreelandella alkaliphila"
    if record.description == "Nil, Halomonas arcis_error": record.description = ">416873:Halomonadaceae:Vreelandella arcis"
    if record.description == "Nil, Halomonas venusta_error": record.description = ">44935:Halomonadaceae:Vreelandella venusta"
    if record.description == "Nil, Halomonas zincidurans_error": record.description = ">1178777:Halomonadaceae:Modicisalibacter zincidurans"
    if record.description == "Nil, Halomonas rituensis_error": record.description = ">2282306:Halomonadaceae:Vreelandella rituensis"
    if record.description == "Nil, Halomonas anticariensis_error": record.description = ">258591:Halomonadaceae:Litchfieldella anticariensis"
    if record.description == "Nil, Halomonas boliviensis_error": record.description = ">223527:Halomonadaceae:Vreelandella boliviensis"
    if record.description == "Nil, Xanthomonas cassavae_error": record.description = ">56450:Xanthomonadaceae:Xanthomonas cassavae"
    if record.description == "Nil, Halomonas subglaciescola_error": record.description = ">29571:Halomonadaceae:Vreelandella subglaciescola"
    if record.description == "Nil, Halomonas endophytica_error": record.description = ">2033802:Halomonadaceae:Billgrantia endophytica"
    if record.description == "Nil, Halomonas desiderata_error": record.description = ">52021:Halomonadaceae:Billgrantia desiderata"
    if record.description == "Nil, Halomonas lutea_error": record.description = ">453962:Halomonadaceae:Modicisalibacter luteus"
    if record.description == "Nil, Halomonas pacifica_error": record.description = ">77098:Halomonadaceae:Bisbaumannia pacifica"
    if record.description == "Nil, Halomonas coralii_error": record.description = ">2304602:Halomonadaceae:Modicisalibacter coralii"
    if record.description == "Nil, Halomonas qijiaojingensis_error": record.description = ">980347:Halomonadaceae:Litchfieldella qijiaojingensis"
    if record.description == "Nil, Halomonas pellis_error": record.description = ">2606936:Halomonadaceae:Billgrantia pellis"
    if record.description == "Nil, Halomonas xianhensis_error": record.description = ">442341:Halomonadaceae:Modicisalibacter xianhensis"
    if record.description == "Nil, Halomonas salicampi_error": record.description = ">1449798:Halomonadaceae:Vreelandella salicampi"
    if record.description == "Nil, Halomonas azerica_error": record.description = ">2732867:Halomonadaceae:Vreelandella azerica"
    if record.description == "Nil, Halomonas saliphila_error": record.description = ">1848458:Halomonadaceae:Billgrantia saliphila"
    if record.description == "Nil, Halomonas lactosivorans_error": record.description = ">2185141:Halomonadaceae:Billgrantia lactosivorans"
    if record.description == "Nil, Halomonas xinjiangensis_error": record.description = ">1166948:Halomonadaceae:Litchfieldella xinjiangensis"
    if record.description == "Nil, Halomonas ilicicola_error": record.description = ">480814:Halomonadaceae:Modicisalibacter ilicicola"
    if record.description == "Nil, Pseudomonas cremoricolorata_error": record.description = ">157783:Pseudomonadaceae:Pseudomonas cremoricolorata"
    if record.description == "Nil, Luteimonas arsenica_error": record.description = ">1586242:Xanthomonadaceae:Luteimonas arsenica"
    if record.description == "Nil, Actinomadura parmotrematis_error": record.description = ">2864039:Thermomonosporaceae:Actinomadura parmotrematis"
    if record.description == "Nil, Neolewinella lacunae_error": record.description = ">1517758:Lewinellaceae:Neolewinella lacunae"
    if record.description == "Nil, Gracilibacillus oryzae_error": record.description = ">1672701:Bacillaceae:Gracilibacillus oryzae"
    if record.description == "Nil, Tetrasphaera jenkinsii_error": record.description = ">330834:Intrasporangiaceae:Nostocoides jenkinsii"
    if record.description == "Nil, Tetrasphaera australiensis_error": record.description = ">99480:Intrasporangiaceae:Nostocoides australiense"
    if record.description == "Nil, Quadrisphaera granulorum_error": record.description = ">317664:Kineosporiaceae:Quadrisphaera granulorum"
    if record.description == "Nil, Rathayibacter oskolensis_error": record.description = ">1891671:Microbacteriaceae:Rathayibacter oskolensis"
    if record.description == "Nil, Paenibacillus vietnamensis_error": record.description = ">2590547:Paenibacillaceae:Paenibacillus vietnamensis"
    if record.description == "Nil, Lawsonibacter celer_error": record.description = ">2986526:Oscillospiraceae:Lawsonibacter celer"
    if record.description == "Nil, Phormidium tenue_error": record.description = ">126344:Oscillatoriaceae:Phormidium tenue"
    if record.description == "Nil, Virgibacillus halodenitrificans_error": record.description = "1482:Bacillaceae:Virgibacillus halodenitrificans"

# Filtering sequences not appropriate to keep in sequence set
count_removed = 0
count_keep = 0

count_nil=0
count_envi=0
count_uncul=0
count_unclass=0
count_candidate=0
count_partial=0

for record in records:
    # Removed: Sequences which failed API calls to Taxallnomy and could not be manually replaced as above
    if "Nil," in record.description: # 16 sequences removed
        remove.append(record)
        count_removed+=1
        count_nil+=1
    # Removed: Sequences with incomplete taxonomic ID/lineage
    elif "Fam_of_environmental samples" in record.description:
        remove.append(record)
        count_removed+=1
        count_envi+=1
    elif "unclassified" in record.description:
        remove.append(record)
        count_removed+=1
        count_unclass+=1
    # Removed: candidate and uncultured sequences as cannot confirm species with certainty
    elif "andidat" in record.description:
        remove.append(record)
        count_removed+=1
        count_candidate+=1
    elif "uncultured" in record.description:
        remove.append(record)
        count_removed+=1
        count_uncul+=1
    # Removed: partial sequences (-3 S.D. from sample mean based on amidase_3 length of the 19 candidate sequences)
    elif len(record.seq) <= 131:
        remove.append(record)
        count_removed+=1
        count_partial+=1
    # Keep all other sequences which pass filtering checks
    else:
        keep.append(record)
        count_keep+=1

print(count_nil)
print(count_envi)
print(count_uncul)
print(count_unclass)
print(count_candidate)
print(count_partial)

print(count_removed)
print(count_keep)

SeqIO.write(keep, "./01_DATA/Amidase_3/03_2_Sequence_Annotation/keep_seq_taxid.txt", "fasta")
SeqIO.write(remove, "./01_DATA/Amidase_3/03_2_Sequence_Annotation/removed_seq_taxid.txt", "fasta")