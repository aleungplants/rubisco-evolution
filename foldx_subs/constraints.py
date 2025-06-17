input_file = "constraintsdraft.txt"
output_file = "constraints.txt"

chains = ["A","B","C","D","E","F","G","H"]

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:

        if "#" in line:
            continue
        else:
            for chain in chains:

                site = line.strip()
                res = site[0]
                pos = site[1:4]

                outline = res + chain + pos + ";\n"

                outfile.write(outline)