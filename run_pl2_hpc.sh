#!/bin/bash
#SBATCH --job-name=metagenomics_pl2
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

set -eo pipefail

echo "Running on $(hostname) at $(date)"

# ---- mamba/conda initialization ----
eval "$(mamba shell hook --shell bash)"
mamba activate genomad

# ---- paths ----
BASE=/mnt/hpc_acegid/nfsscratch/soumareh
OUTPUT=${BASE}/kraken_summary_files   # kraken reports live here; outputs go here too
FASTQ_DIR=${BASE}/fastq               # directory with *_R1*.fastq.gz files (for sample discovery)

# ---- run pipeline ----
run_metagenomics_pl2 \
  --input_dir   "${FASTQ_DIR}" \
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
  --parallel    8 \
  --max_assemblies 2 \
  --bwa_threads 4 \
  --max_workers 8

echo "Pipeline finished at $(date)"
