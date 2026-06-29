#!/bin/bash
#SBATCH --job-name=metag_phase2
#SBATCH --cpus-per-task=32
#SBATCH --mem=496G
#SBATCH --time=48:00:00
#SBATCH --output=logs/phase2_%j.out
#SBATCH --error=logs/phase2_%j.err

set -eo pipefail

echo "Phase 2 on $(hostname) at $(date)"

# ---- mamba/conda ----
eval "$(mamba shell hook --shell bash)"
mamba activate genomad

# ---- paths ----
RUNDIR=/mnt/hpc_acegid/home/soumareh/260424_VH00635_6_AAF2HNTHV
FASTQ_DIR=${RUNDIR}/fastq
BASE=/mnt/hpc_acegid/nfsscratch/soumareh
OUTPUT_DIR=${BASE}/kraken_summary_files

# ---- collect all per-sample contigs into shared output dir ----
# The pipeline expects {sample}_contigs.fasta flat in output_dir
echo "Collecting contigs from per-sample directories..."
for SAMPLE_OUT in "${OUTPUT_DIR}"/*/; do
    SAMPLE=$(basename "${SAMPLE_OUT}")
    SRC="${SAMPLE_OUT}/${SAMPLE}_contigs.fasta"
    DST="${OUTPUT_DIR}/${SAMPLE}_contigs.fasta"
    if [[ -f "${SRC}" && ! -f "${DST}" ]]; then
        cp "${SRC}" "${DST}"
        echo "  Collected: ${SAMPLE}_contigs.fasta"
    fi
done

# ---- collect all Kraken reports into shared output dir ----
for SAMPLE_OUT in "${OUTPUT_DIR}"/*/; do
    SAMPLE=$(basename "${SAMPLE_OUT}")
    SRC="${SAMPLE_OUT}/${SAMPLE}_kraken_report.txt"
    DST="${OUTPUT_DIR}/${SAMPLE}_kraken_report.txt"
    if [[ -f "${SRC}" && ! -f "${DST}" ]]; then
        cp "${SRC}" "${DST}"
    fi
done

echo "Contigs collected:"
ls "${OUTPUT_DIR}"/*_contigs.fasta | wc -l

# ---- Phase 2: geNomad + DIAMOND + alignment on all samples ----
run_metagenomics_pl2 \
  --input_dir        "${FASTQ_DIR}" \
  --output_dir       "${OUTPUT_DIR}" \
  --no_metadata \
  --threads          "${SLURM_CPUS_PER_TASK}" \
  --kraken_db        /mnt/hpc_acegid/nfsscratch/DATABASE/Kraken \
  --use_precomputed_reports \
  --skip_existing \
  --no_bowtie2 \
  --diamond \
  --diamond_db       "${BASE}/MNT/nr_genomad.dmnd" \
  --genomad_db       "${BASE}/genomad_db/" \
  --nr_path          "${BASE}/MNT/nr" \
  --run_alignment \
  --skip_reports \
  --skip_multiqc \
  --parallel         8 \
  --max_assemblies   1 \
  --bwa_threads      4 \
  --max_workers      8

echo "Phase 2 done at $(date)"
