#!/bin/bash

# cd the folder containing this file
cd "$(dirname "$0")"

mkdir -p site
cd site

while IFS=, read -r line; do
  csubst site \
  --alignment_file ../clean_aligned_csubst.fasta \
  --rooted_tree_file ../dryopteris_csubst.tre \
  --cb_file ../branch/csubst_cb_2.tsv \
  --branch_id $line \
  --genetic_code 11 \
  --threads 8
done < ../branch/branches_of_interest.txt