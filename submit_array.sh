#!/bin/bash
# Run this script ONCE to generate the sample list and submit the array job.
# Usage: bash submit_array.sh

BASE=/mnt/hpc_acegid/nfsscratch/soumareh
FASTQ_ROOT=${BASE}/fastq          # root dir — one subdirectory per sample

# ---- generate sample list (one sample dir per line) ----
find "${FASTQ_ROOT}" -mindepth 1 -maxdepth 1 -type d | sort > "${BASE}/sample_list.txt"

N=$(wc -l < "${BASE}/sample_list.txt")
echo "Found ${N} samples:"
cat "${BASE}/sample_list.txt"

if [[ $N -eq 0 ]]; then
    echo "No samples found. Check FASTQ_ROOT path."
    exit 1
fi

mkdir -p "${BASE}/logs"

# ---- submit array (0-indexed, so last index = N-1) ----
sbatch --array=0-$((N - 1))%8 run_pl2_array.sh
