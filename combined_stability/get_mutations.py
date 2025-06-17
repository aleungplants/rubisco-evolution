from Bio import AlignIO
from Bio import PDB
import pandas as pd

# Get the aligned gene sequences
with open("complete_aligned_spinach_checked.fas", "r") as infile:
    data = AlignIO.read(infile, "fasta")

aa = []
for i in range(len(data)):
    species_aa = [data[i].name]
    sequence = data[i].seq.translate(gap = "-", table = 11)
    species_aa.append(str(sequence))
    aa.append(species_aa)

print("Translated gene sequences")

# Ensure unique values while preserving order
seqs_only = list(dict.fromkeys([row[1] for row in aa]))
print(len(seqs_only), "unique amino acid sequences")

# Write to a file
with open("seqs_in_order.txt", "w") as f:
    f.writelines(seq + "\n" for seq in seqs_only)

# Get the protein sequence from the PDB
pdb_file = "../foldx_subs/8ruc-assembly-clean_Repair.pdb"
parser = PDB.PDBParser(QUIET = True)
structure = parser.get_structure("protein", pdb_file)

three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

seq_8RUC_res = []
for model in structure:
    for chain in model:
        if chain.id == "A":  # Select chain A
            print(chain)
            for residue in chain:
                res_name = residue.get_resname()
                if res_name in three_to_one:  # Ignore non-standard residues
                    seq_8RUC_res.append(three_to_one[res_name])
                else:
                    print(res_name, "is a nonstandard amino acid")

# Convert list to string
seq_8RUC = "".join(seq_8RUC_res)
print(seq_8RUC)

print("Got PDB protein sequence")

seqs_fixed = []

for seq in seqs_only:
    seq = seq[8:] # 8RUC only has site 9 onwards
    print(seq.count("-"), "ambiguous sites")
    seq_fixed = ""
    for j in range(len(seq)):
        if seq[j] == "-":
            seq_fixed = seq_fixed + seq_8RUC[j]
        else:
            seq_fixed = seq_fixed + seq[j]
    seqs_fixed.append(seq_fixed)

with open("mutant_file.txt", "w") as f:
    for seq in seqs_fixed:
        f.write(seq + "\n")

chains = ["A", "B", "C", "D", "E", "F", "G", "H"]
individual_list = []

count = 0

for seq in seqs_fixed:
    seq_subs = []
    for k in range(len(seq)):
        if seq[k] != seq_8RUC[k]:
            for chain in chains:
                sub = seq_8RUC[k] + chain + str(k + 9) + seq[k]
                seq_subs.append(sub)
                count = count + 1
    seq_subs = ",".join(seq_subs) + ";"
    individual_list.append(seq_subs)

print("Total number of mutations:", count)

with open("individual_list_mutations.txt", "w") as f:
    for mutations in individual_list:
        f.write(mutations + "\n")

# Split into 8 chunks
input_file = "individual_list_mutations.txt"
with open(input_file, "r") as file:
    lines = file.readlines()

total_lines = len(lines) # Calculate the number of lines per chunk and the remainder (if any)
print("Read", total_lines, "lines")

num_chunks = 8
chunk_size = total_lines // num_chunks
remainder = total_lines % num_chunks

chunks = [] # Split lines into chunks
start_index = 0
for i in range(num_chunks):
    # Distribute the remainder by adding 1 line to the first 'remainder' chunks
    end_index = start_index + chunk_size + (1 if i < remainder else 0)
    chunks.append(lines[start_index:end_index])
    print("Got chunk", str(i + 1), "with lines", start_index, "to", str(end_index - 1), ", inclusive")
    start_index = end_index

# Write each chunk to a new file
for i, chunk in enumerate(chunks, 1):
    output_filename = f"{input_file.rstrip('.txt')}_{i}.txt"
    with open(output_filename, "w") as output_file:
        output_file.writelines(chunk)
    print(f"Created: {output_filename}")