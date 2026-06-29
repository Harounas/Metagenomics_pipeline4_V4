#!/bin/bash
#SBATCH --job-name=viral_classification
#SBATCH --cpus-per-task=32
#SBATCH --mem=496G
#SBATCH --time=72:00:00
#SBATCH --output=logs/phase3_%j.out
#SBATCH --error=logs/phase3_%j.err

# Phase 3: Viral classification from merged contigs
#
# Input:  per-sample *_contigs.fasta files already present in OUTPUT_DIR
#         (produced by phase1_assembly_array.sh)
#
# Steps performed by run_viral_classification:
#   1. Merge all per-sample contigs -> merged_contigs.fasta
#   2. Kraken2 on merged contigs    -> kraken_output/
#   3. DIAMOND (initial pass)       -> diamond_initial.m8
#   4. Extract viral contigs        -> viral_contigs_merged.fasta
#   5. CD-HIT-EST clustering        -> clustered_output/
#   6. DIAMOND on clustered contigs -> results_clustered.m8
#   7. Virus-name annotation        -> diamond_results_contig_with_sampleid.tsv
#   8. Final TSV                    -> filtered_clusters_assigned_rep_virus.tsv
#
# Usage:
#   sbatch phase3_viral_classification.sh
#   # or, to chain after phase1:
#   sbatch --dependency=afterok:<PHASE1_JOBID> phase3_viral_classification.sh

set -eo pipefail

echo "Phase 3 starting on $(hostname) at $(date)"

# ---- mamba/conda ----
eval "$(mamba shell hook --shell bash)"
mamba activate genomad

# ---- paths ---- (edit to match your environment) ----
BASE=/mnt/hpc_acegid/nfsscratch/soumareh
OUTPUT_DIR=${BASE}/kraken_summary_files   # must contain *_contigs.fasta files

KRAKEN_DB=/mnt/hpc_acegid/nfsscratch/DATABASE/Kraken
DIAMOND_DB=${BASE}/MNT/nr_genomad.dmnd
NR_PATH=${BASE}/MNT/nr

mkdir -p "${BASE}/logs"

# ---- run viral classification workflow ----
run_viral_classification \
    --output_dir   "${OUTPUT_DIR}" \
    --kraken_db    "${KRAKEN_DB}" \
    --diamond_db   "${DIAMOND_DB}" \
    --nr_path      "${NR_PATH}" \
    --threads      "${SLURM_CPUS_PER_TASK}" \
    --min_length   200 \
    --skip_existing

# ---- results summary ----
FINAL_TSV="${OUTPUT_DIR}/filtered_clusters_assigned_rep_virus.tsv"
if [[ -f "${FINAL_TSV}" ]]; then
    N=$(tail -n +2 "${FINAL_TSV}" | wc -l)
    echo "Phase 3 complete at $(date)"
    echo "Final TSV: ${FINAL_TSV}  (${N} viral contig rows)"
else
    echo "ERROR: Final TSV not produced." >&2
    exit 1
fi
