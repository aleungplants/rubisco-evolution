from Bio import AlignIO

with open("clean_aligned_checked.fas", "r") as infile:
    data = AlignIO.read(infile, "fasta")
    print("Read file with", len(data), "records")

    for item in data:
        name = item.description.split(" ")[1]
        if len(item.description) > 10:
            name = name[0:10]
        if len(item.description.split(" ")) > 2:
            if len(item.description) > 7:
                name = name[0:7]
            name = name + item.description.split(" ")[2][1:4]
        item.id = name
        item.name = name
        # print(len(item._seq))
        item.description = ""





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
