#!/bin/bash

# cd the folder containing this file
cd "$(dirname "$0")"

mkdir -p branch
cd branch

csubst analyze \
--alignment_file ../clean_aligned_csubst.fasta \
--rooted_tree_file ../Pinus_time_csubst.tre \
--genetic_code 11 \
--iqtree_redo yes \
--threads 3


