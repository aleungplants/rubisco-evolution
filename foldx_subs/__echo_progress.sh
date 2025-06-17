#!/bin/bash

njobs_parallel=8

# Ask for start time input if not already supplied
read -p "Enter start time (YYYY-MM-DD HH:MM:SS or 'now'): " start_time

# If input is "now" or empty, use current time
if [[ -z "$start_time" || "$start_time" == "now" ]]; then
  start_epoch=$(date +%s)
else
  # Try converting input to epoch
  if ! start_epoch=$(date -d "$start_time" +%s 2>/dev/null); then
    echo "❌ Invalid date format. Use 'YYYY-MM-DD HH:MM:SS' or type 'now'."
    exit 1
  fi
fi

# Change to the directory where this script is located
cd "$(dirname "$0")"
cd /runs

# Run indefinitely
while true; do
  clear
  now_epoch=$(date +%s)
  elapsed=$((now_epoch - start_epoch))
  elapsed_fmt=$(printf '%02d:%02d:%02d\n' $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60)))

  echo -e "\e[1;34m[$(date '+%Y-%m-%d %H:%M:%S')] FoldX Progress\e[0m"
  echo "Elapsed time: $elapsed_fmt"

  cd runs || exit 1

  mutation_file="../individual_list_mutations.txt"

  if [[ -f "$mutation_file" ]]; then
    TOTAL_ROWS=$(wc -l < "$mutation_file")
  else
    echo "ERROR: Mutation file not found: $mutation_file"
    TOTAL_ROWS=-1
  fi

  SUM_ROWS=0

  for i in $(seq 1 $njobs_parallel); do
    fxout=$(ls Average_run_${i}_*.fxout 2>/dev/null)
    if [[ -n "$fxout" ]]; then
      ROWS=$(tail -n +10 "$fxout" | wc -l)
      SUM_ROWS=$((SUM_ROWS + ROWS))
      echo "run_$i -> $ROWS rows"
    else
      echo "run_$i -> (no file yet)"
    fi
  done

  if [[ $TOTAL_ROWS -gt 0 ]]; then
    PERCENT=$(echo "scale=3; 100 * $SUM_ROWS / $TOTAL_ROWS" | bc)
    INT_PERCENT=$(echo "$PERCENT / 1" | bc)
  else
    PERCENT="0.000"
    INT_PERCENT=0
  fi

  # Draw progress bar
  BAR_WIDTH=30
  FILLED=$(( (INT_PERCENT * BAR_WIDTH) / 100 ))
  EMPTY=$(( BAR_WIDTH - FILLED ))

  printf "Progress: ["
  printf "%0.s█" $(seq 1 $FILLED)
  printf "%0.s░" $(seq 1 $EMPTY)
  printf "] %.3f%% (%d / %d)\n" "$PERCENT" "$SUM_ROWS" "$TOTAL_ROWS"

  # Estimate ETA
  if (( INT_PERCENT > 0 )); then
    est_total_sec=$(echo "$elapsed * 100 / $PERCENT" | bc -l | cut -d'.' -f1)
    eta_sec=$((est_total_sec - elapsed))

    days=$((eta_sec / 86400))
    hours=$(( (eta_sec % 86400) / 3600 ))
    minutes=$(( (eta_sec % 3600) / 60 ))
    seconds=$((eta_sec % 60))

    if (( days > 0 )); then
      eta_fmt=$(printf '%dd %02d:%02d:%02d\n' "$days" "$hours" "$minutes" "$seconds")
    else
      eta_fmt=$(printf '%02d:%02d:%02d\n' "$hours" "$minutes" "$seconds")
    fi

    echo "ETA (estimated remaining time): $eta_fmt"
  else
    echo "ETA: --:--:-- (waiting for progress)"
  fi

  COUNT=$(pgrep -a foldx | grep -v 'grep' | wc -l)
  echo "Running FoldX processes: $COUNT / $njobs_parallel"

  cd .. || exit 1
  sleep 600
done
