import csv
from Bio import Entrez, SeqIO

sessa_species = [] # Sessa et al. 2017 AJB
# File containing accessions from Sessa et al. 2017 AJB and species names. Names are needed because GenBank names are not always accurate.
reader = csv.DictReader(open("sessa_tree_accessions.csv", encoding = "utf-8-sig"))
accessions = {}
for row in reader:
    accessions[row["Taxon"]] = row["rbcL"]
    sessa_species.append(row["Taxon"])

# Download from the list of accessions
GB_FILENAME = "sequence.gb"
Entrez.email = "arthur.aqualung@gmail.com"
Entrez.api_key = "74e17638c83e6acc5985a5948e22c318a708"

open(GB_FILENAME, "w")
for species, accession in accessions.items():
    if accession == "—":
        continue
    record = Entrez.efetch(db = "nucleotide", id = accession, rettype = "gb", retmode = "text")
    print("Writing", accession, "to", GB_FILENAME)
    with open(GB_FILENAME, "a") as f:
        f.write(record.read())
print("Wrote .gb file")
