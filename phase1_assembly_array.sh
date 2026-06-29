#!/bin/bash
#SBATCH --job-name=metag_phase1
#SBATCH --cpus-per-task=32
#SBATCH --mem=496G
#SBATCH --time=48:00:00
#SBATCH --array=0-999%8          # adjust %8 = max simultaneous jobs; 0-999 = max samples
#SBATCH --output=logs/phase1_%A_%a.out
#SBATCH --error=logs/phase1_%A_%a.err

set -eo pipefail

echo "Task ${SLURM_ARRAY_TASK_ID} on $(hostname) at $(date)"

# ---- mamba/conda ----
eval "$(mamba shell hook --shell bash)"
mamba activate genomad

# ---- paths ----
# FASTQ source: login-node path (may contain symlinks to acegid_storage)
RUNDIR=/mnt/hpc_acegid/home/soumareh/260424_VH00635_6_AAF2HNTHV
FASTQ_SRC=${RUNDIR}/fastq          # accessible on login node only
BASE=/mnt/hpc_acegid/nfsscratch/soumareh
OUTPUT_DIR=${BASE}/kraken_summary_files

# ---- resolve actual FASTQ path (follow symlinks on login node at submit time) ----
# If files were already copied to nfsscratch use that, else abort with clear message
FASTQ_DIR=${BASE}/fastq
if [[ ! -d "${FASTQ_DIR}" ]]; then
    echo "ERROR: ${FASTQ_DIR} not accessible on $(hostname)"
    echo "Run on login node first:"
    echo "  rsync -avL ${FASTQ_SRC}/ ${FASTQ_DIR}/"
    exit 1
fi

# ---- get this task's sample directly from FASTQ dir (no sample_list.txt needed) ----
SAMPLE=$(ls "${FASTQ_DIR}"/*_R1*.fastq.gz 2>/dev/null \
         | sort \
         | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" \
         | xargs basename 2>/dev/null \
         | sed 's/_R1_001\.fastq\.gz$//' \
         | sed 's/_R1\.fastq\.gz$//')

if [[ -z "${SAMPLE}" ]]; then
    echo "No sample for task ${SLURM_ARRAY_TASK_ID} — exiting."
    exit 0
fi
echo "Sample: ${SAMPLE}"

# ---- create per-sample input dir with symlinks ----
SAMPLE_INPUT=${SLURM_TMPDIR}/${SAMPLE}_input
mkdir -p "${SAMPLE_INPUT}"

# Symlink R1 and R2 for this sample (handles _R1.fastq.gz and _R1_001.fastq.gz)
for R1 in "${FASTQ_DIR}/${SAMPLE}"*_R1*.fastq.gz; do
    [[ -e "$R1" ]] && ln -sf "$R1" "${SAMPLE_INPUT}/"
done
for R2 in "${FASTQ_DIR}/${SAMPLE}"*_R2*.fastq.gz; do
    [[ -e "$R2" ]] && ln -sf "$R2" "${SAMPLE_INPUT}/"
done

echo "Input files:"
ls -lh "${SAMPLE_INPUT}/"

# ---- create per-sample output dir ----
SAMPLE_OUT=${OUTPUT_DIR}/${SAMPLE}
mkdir -p "${SAMPLE_OUT}"

# ---- copy precomputed Kraken report if it exists ----
KREPORT="${BASE}/kraken_summary_files/${SAMPLE}_kraken_report.txt"
if [[ -f "${KREPORT}" && ! -f "${SAMPLE_OUT}/${SAMPLE}_kraken_report.txt" ]]; then
    cp "${KREPORT}" "${SAMPLE_OUT}/${SAMPLE}_kraken_report.txt"
fi

# ---- copy existing contigs if they exist ----
for CONTIGS_SRC in \
    "${BASE}/kraken_summary_files/${SAMPLE}_contigs.fasta" \
    "${RUNDIR}/${SAMPLE}/contigs.fasta"; do
    CONTIGS_DST="${SAMPLE_OUT}/${SAMPLE}_contigs.fasta"
    if [[ -f "${CONTIGS_SRC}" && ! -f "${CONTIGS_DST}" ]]; then
        cp "${CONTIGS_SRC}" "${CONTIGS_DST}"
        echo "Copied contigs: ${CONTIGS_SRC} → ${CONTIGS_DST}"
        break
    fi
done

# ---- Phase 1: fastp + bowtie2 + SPAdes + Kraken (no DIAMOND/geNomad yet) ----
run_metagenomics_pl2 \
  --input_dir        "${SAMPLE_INPUT}" \
  --output_dir       "${SAMPLE_OUT}" \
  --no_metadata \
  --threads          "${SLURM_CPUS_PER_TASK}" \
  --kraken_db        /mnt/hpc_acegid/nfsscratch/DATABASE/Kraken \
  --bowtie2_index    /mnt/hpc_acegid/home/soumareh/haouruna/GR38_bt2 \
  --use_precomputed_reports \
  --use_assembly \
  --skip_existing \
  --skip_reports \
  --skip_multiqc \
  --parallel         1 \
  --max_assemblies   1

echo "Phase 1 done for ${SAMPLE} at $(date)"
