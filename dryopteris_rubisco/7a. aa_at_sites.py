import pandas as pd
from Bio import AlignIO


with open("clean_aligned_csubst.fasta", "r") as infile:
    data = AlignIO.read(infile, "fasta")

omega = pd.read_csv("pos_selected_conv_sites.csv")
sites = omega["Site"].unique()
sites = list(sites)
print(sites)

aa = []
for i in range(len(data)):
    species_aa = [data[i].name]
    for j in range(len(sites)):
        print(sites[j])
        start = 3 * (sites[j] - 1 - 20)
        end = 3 * (sites[j] - 1 - 20) + 3
        try:
            species_aa.append(data[i].seq[start:end].translate(gap = "-", table = 11))
            print(species_aa)
        except:
            print("Error with translating at taxon", i,"and AA position", sites[j],
                  ", in", data[i].name, ".",
                  "The codon was", data[i].seq[start:end])
            species_aa.append("-")
            continue
    aa.append(species_aa)

headers = sites
headers.insert(0, "Species")
aa = pd.DataFrame(aa, columns = headers)
aa.to_csv("aa_at_sites.csv", index = False)
