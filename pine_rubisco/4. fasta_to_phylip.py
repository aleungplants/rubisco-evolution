from Bio import AlignIO

with open("aligned_checked.fasta", "r") as infile:
    data = AlignIO.read(infile, "fasta")
    print("Read file with", len(data), "records")

metadata = ["Species,Accession\n"]
list = []
for item in data:
    species = " ".join(item.description.split(" ")[0:2]) # description is the fasta entry name, Genus epithet ACCESSION
    if len(item.description.split(" ")) > 3:
        accession = "_".join(item.description.split(" ")[2:])
    else:
        accession = item.description.split(" ")[2]
    metadata.append(species + "," + accession + "\n")
    item.description = item.description.split(" ")[1]
    if len(item.description) > 10:
        item.description = item.description[0:10]
    item.id = item.description
    item.name = item.description

with open("accessions.csv", 'w') as outfile:
    outfile.writelines(metadata)

with open("clean_aligned.phylip", "w") as outfile:
    AlignIO.write(data, outfile, "phylip")
    print("Saved", "clean_aligned.phylip", "with", len(data), "records")

with open("clean_aligned_csubst.fasta", "w") as outfile:
    AlignIO.write(data, outfile, "fasta")
    print("Saved", "clean_aligned_csubst.fasta", "with", len(data), "records")

# Fix formatting for PAML: add "I" to line 1 and add spaces after taxon names
with open("clean_aligned.phylip", "r") as infile, open("clean_aligned_paml.phylip", "w") as outfile:
    for line_number, line in enumerate(infile):
        if line_number == 0:
            line = line.replace("\n", ' I\n')
        else:
            line = line.split(" ", maxsplit = 1) # split removes the delimiter
            line = "  ".join(line) # add back delimited (space) with an extra space - PAML wants two spaces between taxon name and sequence
        outfile.write(line)
print("Saved", "clean_aligned_paml.phylip", "in PAML-friendly format")
