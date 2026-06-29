#!/bin/bash
# Phase 1: fastp + bowtie2 (host depletion) + SPAdes assembly + Kraken2
# Processes all samples in FASTQ_DIR sequentially.
# Usage: bash run_phase1.sh

set -eo pipefail

# ============================================================
# EDIT THESE PATHS
# ============================================================
BASE=/mnt/hpc_acegid/nfsscratch/soumareh
FASTQ_DIR=${BASE}/fastq
OUTPUT_DIR=${BASE}/kraken_summary_files
KRAKEN_DB=/mnt/hpc_acegid/nfsscratch/DATABASE/Kraken
BOWTIE2_INDEX=/mnt/hpc_acegid/home/soumareh/haouruna/GR38_bt2
THREADS=32
# ============================================================

mkdir -p "${OUTPUT_DIR}" logs

echo "=========================================="
echo "Phase 1 started at $(date)"
echo "FASTQ dir : ${FASTQ_DIR}"
echo "Output dir: ${OUTPUT_DIR}"
echo "=========================================="

# Collect all samples from R1 files
SAMPLES=()
for R1 in "${FASTQ_DIR}"/*_R1*.fastq.gz; do
    [[ -e "$R1" ]] || { echo "ERROR: No *_R1*.fastq.gz files found in ${FASTQ_DIR}"; exit 1; }
    SAMPLE=$(basename "$R1" | sed 's/_R1_001\.fastq\.gz$//' | sed 's/_R1\.fastq\.gz$//')
    SAMPLES+=("$SAMPLE")
done

echo "Found ${#SAMPLES[@]} samples"

for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "------------------------------------------"
    echo "Processing: ${SAMPLE}  ($(date))"
    echo "------------------------------------------"

    SAMPLE_INPUT=$(mktemp -d)
    SAMPLE_OUT=${OUTPUT_DIR}/${SAMPLE}
    mkdir -p "${SAMPLE_OUT}"

    # Symlink FASTQ files into temp input dir
    for R1 in "${FASTQ_DIR}/${SAMPLE}"*_R1*.fastq.gz; do
        [[ -e "$R1" ]] && ln -sf "$R1" "${SAMPLE_INPUT}/"
    done
    for R2 in "${FASTQ_DIR}/${SAMPLE}"*_R2*.fastq.gz; do
        [[ -e "$R2" ]] && ln -sf "$R2" "${SAMPLE_INPUT}/"
    done

    # Copy precomputed Kraken report if it exists
    KREPORT="${BASE}/kraken_summary_files/${SAMPLE}_kraken_report.txt"
    [[ -f "$KREPORT" && ! -f "${SAMPLE_OUT}/${SAMPLE}_kraken_report.txt" ]] && \
        cp "$KREPORT" "${SAMPLE_OUT}/${SAMPLE}_kraken_report.txt"

    # Copy existing contigs if available
    for CONTIGS_SRC in \
        "${BASE}/kraken_summary_files/${SAMPLE}_contigs.fasta" \
        "${FASTQ_DIR}/../${SAMPLE}/contigs.fasta"; do
        CONTIGS_DST="${SAMPLE_OUT}/${SAMPLE}_contigs.fasta"
        if [[ -f "$CONTIGS_SRC" && ! -f "$CONTIGS_DST" ]]; then
            cp "$CONTIGS_SRC" "$CONTIGS_DST"
            echo "  Copied existing contigs: $CONTIGS_SRC"
            break
        fi
    done

    run_metagenomics_pl2 \
        --input_dir              "${SAMPLE_INPUT}" \
        --output_dir             "${SAMPLE_OUT}" \
        --no_metadata \
        --threads                "${THREADS}" \
        --kraken_db              "${KRAKEN_DB}" \
        --bowtie2_index          "${BOWTIE2_INDEX}" \
        --use_precomputed_reports \
        --use_assembly \
        --skip_existing \
        --skip_reports \
        --skip_multiqc \
        --parallel               1 \
        --max_assemblies         1 \
        2>&1 | tee -a "logs/phase1_${SAMPLE}.log"

    rm -rf "${SAMPLE_INPUT}"
    echo "Done: ${SAMPLE}"
done

echo ""
echo "=========================================="
echo "Phase 1 complete at $(date)"
echo "=========================================="
