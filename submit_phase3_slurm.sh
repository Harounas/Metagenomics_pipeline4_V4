#!/bin/bash
#SBATCH --job-name=phase3_viral
#SBATCH --output=/mnt/hpc_acegid/nfsscratch/soumareh/kraken_summary_files/Contigs/logs/phase3_%j.out
#SBATCH --error=/mnt/hpc_acegid/nfsscratch/soumareh/kraken_summary_files/Contigs/logs/phase3_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --partition=normal

set -eo pipefail

# ── paths (edit these if needed) ─────────────────────────────────────────────
OUTPUT_DIR=/mnt/hpc_acegid/nfsscratch/soumareh/kraken_summary_files/Contigs
KRAKEN_DB=/mnt/hpc_acegid/nfsscratch/DATABASE/Kraken
DIAMOND_DB=/mnt/hpc_acegid/nfsscratch/DATABASE/diamond/nr.dmnd
GENOMAD_DB=/mnt/hpc_acegid/home/soumareh/haouruna/genomad_db
THREADS=32
SPLITS=8
GENOMAD_MIN_LENGTH=500

SCRIPT_DIR=/mnt/hpc_acegid/home/soumareh/Metagenomics_pipeline4_V4

# ── environment ───────────────────────────────────────────────────────────────
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate genomad

export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH}"
mkdir -p "${OUTPUT_DIR}/logs"

echo "Job ID      : ${SLURM_JOB_ID}"
echo "Node        : $(hostname)"
echo "Started at  : $(date)"
echo "Output dir  : ${OUTPUT_DIR}"

python "${SCRIPT_DIR}/Metagenomics_pipeline4_V2/viral_classification_workflow.py" \
    --output_dir          "${OUTPUT_DIR}" \
    --kraken_db           "${KRAKEN_DB}" \
    --diamond_db          "${DIAMOND_DB}" \
    --genomad_db          "${GENOMAD_DB}" \
    --threads             "${THREADS}" \
    --splits              "${SPLITS}" \
    --genomad_min_length  "${GENOMAD_MIN_LENGTH}" \
    --skip_existing \
    2>&1 | tee "${OUTPUT_DIR}/logs/phase3.log"

FINAL_TSV="${OUTPUT_DIR}/filtered_clusters_assigned_rep_virus.tsv"
if [[ -f "${FINAL_TSV}" ]]; then
    N=$(tail -n +2 "${FINAL_TSV}" | wc -l)
    echo ""
    echo "Phase 3 complete at $(date)"
    echo "Final TSV : ${FINAL_TSV}"
    echo "Viral rows: ${N}"
else
    echo "ERROR: Final TSV not found at ${FINAL_TSV}" >&2
    exit 1
fi
