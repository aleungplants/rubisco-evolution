from Bio import AlignIO
import csv

folder_list = ["dryopteris_rubisco",
                "limonium_rubisco",
                "pine_rubisco",
                "csubst"]

complete_seqs = AlignIO.MultipleSeqAlignment(None)
complete_seq_lens = []
total_seq_num = 0

for folder in folder_list:
    with open("../" + folder + "/clean_aligned_csubst.fasta", "r") as infile:
        data = AlignIO.read(infile, "fasta")
        seq_num = len(data)
        total_seq_num = total_seq_num + seq_num
        print("Read file with", seq_num, "records")

        for item in data:
            if "-" not in item.seq:
                print(len(item.seq))
                if len(item.seq[0]) < 1425:
                    seq_length = len(item.seq)
                    item.seq = item.seq + "-"*(1425-seq_length)
                complete_seqs.append(item)
                complete_seq_lens.append(seq_length)

print("There are in total", total_seq_num, "sequences.")
print("Out of those, there are", len(complete_seqs), "sequences that do not have gaps.")
print("On average, they are", sum(complete_seq_lens) / len(complete_seq_lens), "bp long.")

with open("complete_seqs.fasta", "w") as outfile:
    AlignIO.write(complete_seqs, outfile, "fasta")
    print("Saved", "complete_seqs.fasta", "with", len(complete_seqs), "records")

# with open("complete_aligned_spinach_checked.fas", "r") as infile:
#     data = AlignIO.read(infile, "fasta")
#     print("Read file with", len(data), "records")
#
# with open("complete_aligned_spinach_checked.phylip", "w") as outfile:
#     AlignIO.write(data, outfile, "phylip")
#     print("Saved", "complete_aligned_spinach_checked.phylip", "with", len(data), "records")

with open("combined_species_list.csv", mode="r") as file:
    reader = csv.reader(file)
    # Skip the header if present
    next(reader)
    genus_species = {rows[1]: rows[0] for rows in reader}
print(genus_species)

with open("complete_aligned_spinach_checked.fas", "r") as infile:
    data = AlignIO.read(infile, "fasta")
    print("Read file with", len(data), "records")

complete_seqs_genus = AlignIO.MultipleSeqAlignment(None)

for item in data:
    if item.name in genus_species:
        item.id  = genus_species[item.name][0] + "_" + item.name
        item.name = ""
        complete_seqs_genus.append(item)

print(complete_seqs_genus)

with open("complete_seqs_genus.fasta", "w") as outfile:
    AlignIO.write(complete_seqs_genus, outfile, "fasta")
    print("Saved", "complete_seqs_genus.fasta", "with", len(complete_seqs_genus), "records")
