#!/bin/bash
#SBATCH --job-name=metag_pl2
#SBATCH --cpus-per-task=32
#SBATCH --mem=496G               # close to node max (504G), leaves 8G for OS
#SBATCH --time=48:00:00
#SBATCH --array=0-999%8          # up to 1000 samples; max 8 running at once — adjust %N to your cluster
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err

set -eo pipefail

echo "Array task ${SLURM_ARRAY_TASK_ID} on $(hostname) at $(date)"

# ---- mamba/conda initialization ----
eval "$(mamba shell hook --shell bash)"
mamba activate genomad

# ---- paths ----
BASE=/mnt/hpc_acegid/nfsscratch/soumareh
SAMPLE_LIST=${BASE}/sample_list.txt   # one sample directory path per line (generated below)
SAMPLE_DIR=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${SAMPLE_LIST}")

if [[ -z "${SAMPLE_DIR}" ]]; then
    echo "No sample for task ${SLURM_ARRAY_TASK_ID} — exiting."
    exit 0
fi

SAMPLE=$(basename "${SAMPLE_DIR}")
OUTPUT=${BASE}/kraken_summary_files/${SAMPLE}
mkdir -p "${OUTPUT}"

echo "Processing sample: ${SAMPLE}"
echo "Input:  ${SAMPLE_DIR}"
echo "Output: ${OUTPUT}"

# ---- copy existing contigs to output dir if needed ----
# Expected name: {SAMPLE}_contigs.fasta
CONTIGS_SRC="${SAMPLE_DIR}/contigs.fasta"
CONTIGS_DST="${OUTPUT}/${SAMPLE}_contigs.fasta"
if [[ -f "${CONTIGS_SRC}" && ! -f "${CONTIGS_DST}" ]]; then
    cp "${CONTIGS_SRC}" "${CONTIGS_DST}"
    echo "Copied contigs.fasta → ${CONTIGS_DST}"
fi

# Copy precomputed Kraken report if it exists in the source dir
KREPORT_SRC="${SAMPLE_DIR}/${SAMPLE}_kraken_report.txt"
KREPORT_DST="${OUTPUT}/${SAMPLE}_kraken_report.txt"
if [[ -f "${KREPORT_SRC}" && ! -f "${KREPORT_DST}" ]]; then
    cp "${KREPORT_SRC}" "${KREPORT_DST}"
fi

# ---- run pipeline for this single sample ----
run_metagenomics_pl2 \
  --input_dir   "${SAMPLE_DIR}" \
  --output_dir  "${OUTPUT}" \
  --no_metadata \
  --threads     "${SLURM_CPUS_PER_TASK}" \
  --kraken_db   /mnt/hpc_acegid/nfsscratch/DATABASE/Kraken \
  --use_precomputed_reports \
  --skip_existing \
  --no_bowtie2 \
  --diamond \
  --diamond_db  "${BASE}/MNT/nr_genomad.dmnd" \
  --genomad_db  "${BASE}/genomad_db/" \
  --nr_path     "${BASE}/MNT/nr" \
  --run_alignment \
  --parallel    1 \
  --max_assemblies 1 \
  --bwa_threads "${SLURM_CPUS_PER_TASK}" \
  --max_workers 8

echo "Sample ${SAMPLE} done at $(date)"
