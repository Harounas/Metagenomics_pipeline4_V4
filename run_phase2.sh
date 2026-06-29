#!/bin/bash
# Phase 2: collect contigs, run geNomad + DIAMOND + alignment
# Usage: bash run_phase2.sh

set -eo pipefail

# ============================================================
# EDIT THESE PATHS
# ============================================================
BASE=/mnt/hpc_acegid/nfsscratch/soumareh
FASTQ_DIR=${BASE}/fastq
OUTPUT_DIR=${BASE}/kraken_summary_files
KRAKEN_DB=/mnt/hpc_acegid/nfsscratch/DATABASE/Kraken
DIAMOND_DB=${BASE}/MNT/nr_genomad.dmnd
GENOMAD_DB=${BASE}/genomad_db
NR_PATH=${BASE}/MNT/nr
THREADS=32
# ============================================================

mkdir -p logs

echo "=========================================="
echo "Phase 2 started at $(date)"
echo "Output dir: ${OUTPUT_DIR}"
echo "=========================================="

# Collect per-sample contigs into flat output dir
echo "Collecting contigs from per-sample directories..."
for SAMPLE_OUT in "${OUTPUT_DIR}"/*/; do
    SAMPLE=$(basename "${SAMPLE_OUT}")
    SRC="${SAMPLE_OUT}/${SAMPLE}_contigs.fasta"
    DST="${OUTPUT_DIR}/${SAMPLE}_contigs.fasta"
    if [[ -f "$SRC" && ! -f "$DST" ]]; then
        cp "$SRC" "$DST"
        echo "  Collected: ${SAMPLE}_contigs.fasta"
    fi
done

# Collect Kraken reports
for SAMPLE_OUT in "${OUTPUT_DIR}"/*/; do
    SAMPLE=$(basename "${SAMPLE_OUT}")
    SRC="${SAMPLE_OUT}/${SAMPLE}_kraken_report.txt"
    DST="${OUTPUT_DIR}/${SAMPLE}_kraken_report.txt"
    if [[ -f "$SRC" && ! -f "$DST" ]]; then
        cp "$SRC" "$DST"
    fi
done

N=$(ls "${OUTPUT_DIR}"/*_contigs.fasta 2>/dev/null | wc -l)
echo "Total contigs files: ${N}"

run_metagenomics_pl2 \
    --input_dir              "${FASTQ_DIR}" \
    --output_dir             "${OUTPUT_DIR}" \
    --no_metadata \
    --threads                "${THREADS}" \
    --kraken_db              "${KRAKEN_DB}" \
    --use_precomputed_reports \
    --skip_existing \
    --no_bowtie2 \
    --diamond \
    --diamond_db             "${DIAMOND_DB}" \
    --genomad_db             "${GENOMAD_DB}" \
    --nr_path                "${NR_PATH}" \
    --run_alignment \
    --skip_reports \
    --skip_multiqc \
    --parallel               8 \
    --max_assemblies         1 \
    --bwa_threads            4 \
    --max_workers            8 \
    2>&1 | tee logs/phase2.log

echo ""
echo "=========================================="
echo "Phase 2 complete at $(date)"
echo "=========================================="
