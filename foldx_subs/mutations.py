chains = ["A", "B", "C", "D", "E", "F", "G", "H"]
input_file = "individual_listdraft.txt"
output_file = "individual_list_mutations.txt"

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:

        if line[0] == "#" or line == "\n":
            continue

        sub = line.strip()
        pos = sub[1:-1]
        aa = sub[-1]
        aa_anc = sub[0]

        outline = ""
        for i in range(len(chains) - 1):
            outline = outline + aa_anc + chains[i] + pos + aa + ","
        outline = outline + aa_anc + chains[-1] + pos + aa + ";"
        print(outline)
        outfile.write(outline + '\n')


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
