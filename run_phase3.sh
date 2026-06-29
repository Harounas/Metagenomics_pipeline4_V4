#!/bin/bash
# Phase 3: Kraken2 + DIAMOND on merged contigs -> CD-HIT cluster -> final viral TSV
# Usage: bash run_phase3.sh

set -eo pipefail

# ============================================================
# EDIT THESE PATHS
# ============================================================
BASE=/mnt/hpc_acegid/nfsscratch/soumareh
OUTPUT_DIR=${BASE}/kraken_summary_files   # must contain *_contigs.fasta files
KRAKEN_DB=/mnt/hpc_acegid/nfsscratch/DATABASE/Kraken
DIAMOND_DB=${BASE}/MNT/nr_genomad.dmnd
NR_PATH=${BASE}/MNT/nr
THREADS=32
MIN_LENGTH=200
# ============================================================

mkdir -p logs

echo "=========================================="
echo "Phase 3 started at $(date)"
echo "Output dir: ${OUTPUT_DIR}"
echo "=========================================="

run_viral_classification \
    --output_dir     "${OUTPUT_DIR}" \
    --kraken_db      "${KRAKEN_DB}" \
    --diamond_db     "${DIAMOND_DB}" \
    --nr_path        "${NR_PATH}" \
    --threads        "${THREADS}" \
    --min_length     "${MIN_LENGTH}" \
    --skip_existing \
    2>&1 | tee logs/phase3.log

# Check output
FINAL_TSV="${OUTPUT_DIR}/filtered_clusters_assigned_rep_virus.tsv"
if [[ -f "${FINAL_TSV}" ]]; then
    N=$(tail -n +2 "${FINAL_TSV}" | wc -l)
    echo ""
    echo "=========================================="
    echo "Phase 3 complete at $(date)"
    echo "Final TSV : ${FINAL_TSV}"
    echo "Viral rows: ${N}"
    echo "=========================================="
else
    echo "ERROR: Final TSV not found at ${FINAL_TSV}" >&2
    exit 1
fi
