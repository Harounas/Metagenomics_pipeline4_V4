#!/bin/bash
#SBATCH --job-name=metagenomics_e2e
#SBATCH --output=/mnt/hpc_acegid/nfsscratch/soumareh/260424_VH00635_6_AAF2HNTHV/kraken_summary_files/logs/e2e_%j.out
#SBATCH --error=/mnt/hpc_acegid/nfsscratch/soumareh/260424_VH00635_6_AAF2HNTHV/kraken_summary_files/logs/e2e_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=480G
#SBATCH --time=168:00:00
#SBATCH --partition=normal

set -eo pipefail

# ── paths ─────────────────────────────────────────────────────────────────────
FASTQ_DIR=/mnt/hpc_acegid/nfsscratch/soumareh/260424_VH00635_6_AAF2HNTHV/fastq          # <-- update this
OUTPUT_DIR=/mnt/hpc_acegid/nfsscratch/soumareh/260424_VH00635_6_AAF2HNTHV/kraken_summary_files
KRAKEN_DB=/mnt/hpc_acegid/nfsscratch/DATABASE/Kraken
BOWTIE2_INDEX=/mnt/hpc_acegid/home/soumareh/haouruna/GRCh38_bt2   # <-- update this
DIAMOND_DB=/mnt/hpc_acegid/nfsscratch/DATABASE/diamond/nr.dmnd
GENOMAD_DB=/mnt/hpc_acegid/home/soumareh/haouruna/genomad_db
NR_PATH=/mnt/hpc_acegid/nfsscratch/DATABASE/blastdb/nr/nr.faa
SCRIPT_DIR=/mnt/hpc_acegid/home/soumareh/Metagenomics_pipeline4_V4

THREADS=32
MIN_LENGTH=200
GENOMAD_MIN_LENGTH=200
SPLITS=16

# ── environment ───────────────────────────────────────────────────────────────
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate genomad

export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH}"
mkdir -p "${OUTPUT_DIR}/logs"

echo "Job ID      : ${SLURM_JOB_ID}"
echo "Node        : $(hostname)"
echo "Started at  : $(date)"
echo "FASTQ dir   : ${FASTQ_DIR}"
echo "Output dir  : ${OUTPUT_DIR}"

bash "${SCRIPT_DIR}/run_pipeline.sh" \
    --fastq_dir      "${FASTQ_DIR}" \
    --output_dir     "${OUTPUT_DIR}" \
    --kraken_db      "${KRAKEN_DB}" \
    --bowtie2_index  "${BOWTIE2_INDEX}" \
    --diamond_db     "${DIAMOND_DB}" \
    --genomad_db     "${GENOMAD_DB}" \
    --nr_path        "${NR_PATH}" \
    --threads        "${THREADS}" \
    --min_length     "${MIN_LENGTH}" \
    --skip_existing

echo ""
echo "Pipeline finished at $(date)"
