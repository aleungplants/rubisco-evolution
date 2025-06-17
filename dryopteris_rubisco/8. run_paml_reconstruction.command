#!/bin/bash

# cd the folder containing this file
cd "$(dirname "$0")"
cd paml_reconstruction

cp ../clean_aligned_paml.phylip clean_aligned_paml.phylip
cp ../dryopteris_paml.tre dryopteris_paml.tre

./codeml


