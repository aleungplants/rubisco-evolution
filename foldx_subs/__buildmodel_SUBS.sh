#!/bin/bash

# ---- Setup ----

cd "$(dirname "$0")"

njobs_parallel=8

mkdir -p logs runs

# ---- Function to run each FoldX job ----
run_foldx_job() {
  TOKEN="a97yjkw7p92nofyhetgizsrqjjnqmu"
  USER="uqjet2t4fnbafi94bzh4wgvzwyv27x"
  TITLE="FoldX analysis progress"
  i="$1"
  ./foldx.exe \
    --command=BuildModel \
    --pdb=8ruc-assembly-clean_Repair.pdb \
    --mutant-file=individual_list_mutations_$1.txt \
    --fix-residues-file=individual_list_constraints.txt \
    --out-pdb=false \
    --output-dir=runs \
    --output-file=run_"$i" \
    --water=-CRYSTAL \
    --vdwDesign=0 \
    --ionStrength=0.15 \
    --numberOfRuns=3 \
    > "logs/foldx_run_${i}.log" 2>&1 && \
  curl -s -F "token=$TOKEN" -F "user=$USER" -F "title=$TITLE" -F "message=Run $i is done." https://api.pushover.net/1/messages.json
}
export -f run_foldx_job

# ---- Run FoldX jobs in parallel ----
parallel --progress --eta run_foldx_job ::: $(seq 1 $njobs_parallel)

# ---- Clean up ----
echo -e "\n✅ All FoldX jobs are done!"
read -p "Press enter to exit."
