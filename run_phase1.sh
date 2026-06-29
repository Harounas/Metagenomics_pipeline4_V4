#!/bin/bash
# Phase 1: fastp + bowtie2 (host depletion) + SPAdes assembly + Kraken2
#
# Usage:
#   bash run_phase1.sh --fastq_dir <path> --output_dir <path> --kraken_db <path> \
#                      --bowtie2_index <path> [--threads 32]

set -eo pipefail

# ── defaults ────────────────────────────────────────────────────────────────
THREADS=32

usage() {
    echo ""
    echo "Usage: bash run_phase1.sh [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  --fastq_dir      DIR   Directory with *_R1*.fastq.gz files"
    echo "  --output_dir     DIR   Output directory"
    echo "  --kraken_db      DIR   Kraken2 database path"
    echo "  --bowtie2_index  PATH  Bowtie2 host index prefix (without .bt2)"
    echo ""
    echo "Optional:"
    echo "  --threads        INT   CPU threads (default: 32)"
    echo "  --help                 Show this help"
    echo ""
    echo "Example:"
    echo "  bash run_phase1.sh \\"
    echo "    --fastq_dir    /data/fastq \\"
    echo "    --output_dir   /data/output \\"
    echo "    --kraken_db    /db/kraken \\"
    echo "    --bowtie2_index /db/GR38_bt2 \\"
    echo "    --threads 32"
    exit 1
}

# ── parse arguments ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --fastq_dir)      FASTQ_DIR="$2";      shift 2 ;;
        --output_dir)     OUTPUT_DIR="$2";     shift 2 ;;
        --kraken_db)      KRAKEN_DB="$2";      shift 2 ;;
        --bowtie2_index)  BOWTIE2_INDEX="$2";  shift 2 ;;
        --threads)        THREADS="$2";        shift 2 ;;
        --help|-h)        usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# ── validate required arguments ──────────────────────────────────────────────
MISSING=""
[[ -z "$FASTQ_DIR"     ]] && MISSING+="  --fastq_dir\n"
[[ -z "$OUTPUT_DIR"    ]] && MISSING+="  --output_dir\n"
[[ -z "$KRAKEN_DB"     ]] && MISSING+="  --kraken_db\n"
[[ -z "$BOWTIE2_INDEX" ]] && MISSING+="  --bowtie2_index\n"

if [[ -n "$MISSING" ]]; then
    echo "ERROR: Missing required arguments:"
    echo -e "$MISSING"
    usage
fi

[[ ! -d "$FASTQ_DIR"  ]] && { echo "ERROR: --fastq_dir not found: $FASTQ_DIR";  exit 1; }
[[ ! -d "$KRAKEN_DB"  ]] && { echo "ERROR: --kraken_db not found: $KRAKEN_DB";  exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH}"
mkdir -p "${OUTPUT_DIR}" "${OUTPUT_DIR}/logs"

echo "=========================================="
echo "Phase 1 started at $(date)"
echo "FASTQ dir   : ${FASTQ_DIR}"
echo "Output dir  : ${OUTPUT_DIR}"
echo "Kraken DB   : ${KRAKEN_DB}"
echo "Bowtie2 idx : ${BOWTIE2_INDEX}"
echo "Threads     : ${THREADS}"
echo "=========================================="

# ── collect samples ───────────────────────────────────────────────────────────
SAMPLES=()
for R1 in "${FASTQ_DIR}"/*_R1*.fastq.gz; do
    [[ -e "$R1" ]] || { echo "ERROR: No *_R1*.fastq.gz files in ${FASTQ_DIR}"; exit 1; }
    SAMPLE=$(basename "$R1" | sed 's/_R1_001\.fastq\.gz$//' | sed 's/_R1\.fastq\.gz$//')
    SAMPLES+=("$SAMPLE")
done
echo "Found ${#SAMPLES[@]} samples"

# ── process each sample ───────────────────────────────────────────────────────
for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "------------------------------------------"
    echo "Processing: ${SAMPLE}  ($(date))"
    echo "------------------------------------------"

    SAMPLE_INPUT=$(mktemp -d)
    SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"
    mkdir -p "${SAMPLE_OUT}"

    for R1 in "${FASTQ_DIR}/${SAMPLE}"*_R1*.fastq.gz; do
        [[ -e "$R1" ]] && ln -sf "$R1" "${SAMPLE_INPUT}/"
    done
    for R2 in "${FASTQ_DIR}/${SAMPLE}"*_R2*.fastq.gz; do
        [[ -e "$R2" ]] && ln -sf "$R2" "${SAMPLE_INPUT}/"
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
        2>&1 | tee -a "${OUTPUT_DIR}/logs/phase1_${SAMPLE}.log"

    rm -rf "${SAMPLE_INPUT}"
    echo "Done: ${SAMPLE}"
done

echo ""
echo "=========================================="
echo "Phase 1 complete at $(date)"
echo "=========================================="
