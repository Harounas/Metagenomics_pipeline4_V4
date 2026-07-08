#!/bin/bash
#SBATCH --job-name=metagenomics_e2e
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=96:00:00
#SBATCH --partition=normal

# Usage: sbatch --export=CONFIG=projects/myproject.conf submit_e2e.sh
# Or:    CONFIG=projects/myproject.conf bash submit_e2e.sh

set -eo pipefail

[[ -z "$CONFIG" ]] && { echo "ERROR: CONFIG not set. Use: sbatch --export=CONFIG=projects/foo.conf submit_e2e.sh"; exit 1; }
[[ ! -f "$CONFIG" ]] && { echo "ERROR: Config file not found: $CONFIG"; exit 1; }
source "$CONFIG"

# ── shared databases (same for all projects) ─────────────────────────────────
KRAKEN_DB=/mnt/hpc_acegid/nfsscratch/DATABASE/Kraken
DIAMOND_DB=/mnt/hpc_acegid/nfsscratch/DATABASE/diamond/nr.dmnd
GENOMAD_DB=/mnt/hpc_acegid/home/soumareh/haouruna/genomad_db
NR_PATH=/mnt/hpc_acegid/nfsscratch/DATABASE/blastdb/nr/nr.faa
SCRIPT_DIR=/mnt/hpc_acegid/home/soumareh/Metagenomics_pipeline4_V4

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate genomad
export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH}"
mkdir -p "${OUTPUT_DIR}/logs"

#SBATCH --output=${OUTPUT_DIR}/logs/e2e_%j.out
#SBATCH --error=${OUTPUT_DIR}/logs/e2e_%j.err

echo "Job ID      : ${SLURM_JOB_ID}"
echo "Node        : $(hostname)"
echo "Started at  : $(date)"
echo "Config      : ${CONFIG}"
echo "Output dir  : ${OUTPUT_DIR}"

bash "${SCRIPT_DIR}/run_pipeline.sh" \
    --fastq_dir      "${FASTQ_DIR}" \
    --output_dir     "${OUTPUT_DIR}" \
    --kraken_db      "${KRAKEN_DB}" \
    --bowtie2_index  "${BOWTIE2_INDEX}" \
    --diamond_db     "${DIAMOND_DB}" \
    --genomad_db     "${GENOMAD_DB}" \
    --nr_path        "${NR_PATH}" \
    --threads        32 \
    --skip_existing

echo "Pipeline finished at $(date)"
