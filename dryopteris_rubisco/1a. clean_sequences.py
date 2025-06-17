import csv
from Bio import SeqIO

# File containing accessions from Sessa et al. 2017 AJB and species names. Names are needed because GenBank names are not always accurate.
reader = csv.DictReader(open("sessa_tree_accessions.csv", encoding = "utf-8-sig"))
accessions = {}
ids = {}
sessa_species = []
for row in reader:
    if "var." in row["Taxon"] or "subsp." in row["Taxon"]:
        subsp = row["Taxon"].split(" ")[2:4]
        subsp = " " + " ".join(subsp)
    else:
        subsp = ""
    sp = row["Taxon"].split(" ")[0:2]
    sp = " ".join(sp) + subsp
    if row["ID"] != "—":
        ids[row["rbcL"]] = row["ID"]
    accessions[row["rbcL"]] = row["Taxon"]
    sessa_species.append(row["Taxon"])

print(ids)

# Read GenBank records as a SeqIO object.
records = list(SeqIO.parse("sequence.gb", "genbank"))

clean_records = []
record_metadata = [] # to keep track of the original .gb data
gene_lengths = {}

print("Opened .gb file with", len(records), "records")
for i in reversed(range(len(records))):
    keep = False
    for j in range(len(records[i].features)):
        if records[i].features[j].type == "CDS" and "gene" in records[i].features[j].qualifiers:
            if records[i].features[j].qualifiers["gene"] == ["rbcL"]: # check all features for rbcL gene
                cds = records[i].features[j].location # get the location data for the rbcL CDS
                cds_seq = cds.extract(records[i].seq)
                records[i].seq = cds_seq
    sp = accessions[records[i].name]
    if "var." in sp or "subsp." in sp:
        subsp = sp.split(" ")[2:4]
        subsp = " " + " ".join(subsp)
    else:
        subsp = ""
    sp = sp.split(" ")[0:2]
    sp = " ".join(sp)
    if sp not in gene_lengths: # keep the most coverage of rbcL
        keep = True
    elif gene_lengths[sp] < len(records[i].seq):
        keep = True
    elif gene_lengths[sp] == len(records[i].seq):
        taxa_overwrite = ["var. fuscoastra", "var. glabra", "var. affinis"] # if same length, keep the type subspecies
        for taxon in taxa_overwrite:
            if taxon in records[i].description:
                keep = True

    if keep == True:
        gene_lengths[sp] = len(records[i].seq)
        records[i].id = sp + subsp
        if records[i].name in ids:
            records[i].id = records[i].id + " " + ids[records[i].name]
        records[i].description = ""
        clean_records.append(records[i])
        record_metadata.append(",".join([records[i].name, records[i].annotations["organism"], records[i].id]))

# check for duplicate names
clean_names = []

for record in clean_records:
    if record.id not in clean_names:
        clean_names.append(record.id)
    else:
        print(record.id)

# Save files.
# with open("clean.fasta", "w") as outfile:
#     SeqIO.write(clean_records, outfile, "fasta")
# print("Saved .fasta file with", len(clean_records), "sequences")
with open("metadata.csv", "w") as outfile:
    outfile.write("accession,gb_species,species\n")
    for k in range(len(record_metadata)):
        outfile.write(record_metadata[k] + "\n")
print("Saved metadata file with", len(record_metadata), "accessions")
