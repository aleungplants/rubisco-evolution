from Bio import SeqIO
import csv

jin_species = []
with open('jin_species.csv') as csvfile:
    for row in csvfile:
        jin_species.append(row.strip())

gb_file = "sequence.gb"

records = list(SeqIO.parse(gb_file, "genbank"))

records_ids = []
records_new = []

for i in range(len(records)):
    if len(records[i].features) > 2:  # 2 because at least one of type="source" and another type="CDS"
        records[i].description = ""
        records[i].id = records[i].annotations["organism"] + " " + records[i].id
        for j in range(len(records[i].features)):  # check all features for rbcL gene
            if len(records[i].seq) < 1200:
                continue
            if "gene" in records[i].features[j].qualifiers and \
                    records[i].features[j].qualifiers["gene"] == ["rbcL"] and \
                    "product" in records[i].features[j].qualifiers and \
                    "carboxylase" in records[i].features[j].qualifiers["product"][0] and \
                    records[i].features[j].type == "CDS":  # just need the CDS
                # print(records[i].features[j])
                cds = records[i].features[j].location # get the location data for the rbcL CDS
                cds_seq = cds.extract(records[i].seq)
                records[i].seq = cds_seq
                if records[i].annotations["organism"] not in records_ids and \
                        records[i].annotations["organism"] in jin_species:
                    records_ids.append(records[i].annotations["organism"])
                    records_new.append(records[i])
            if len(records[i].seq) < 2000:
                continue

# print(records_ids)
print(len(records_ids))
with open("sequences.fasta", "w") as outfile:
    SeqIO.write(records_new, outfile, "fasta")
