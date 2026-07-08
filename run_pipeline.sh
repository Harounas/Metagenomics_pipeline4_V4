#!/bin/bash
# Master script: run all 3 phases in sequence
#
# Usage:
#   bash run_pipeline.sh --fastq_dir <path> --output_dir <path> \
#                        --kraken_db <path> --bowtie2_index <path> \
#                        --diamond_db <path> --genomad_db <path> \
#                        --nr_path <path> [--threads 32] [--min_length 200]

set -eo pipefail

# ── defaults ─────────────────────────────────────────────────────────────────
THREADS=32
MIN_LENGTH=200
SKIP_EXISTING=false

usage() {
    echo ""
    echo "Usage: bash run_pipeline.sh [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  --fastq_dir      DIR   Directory with *_R1*.fastq.gz files"
    echo "  --output_dir     DIR   Output directory"
    echo "  --kraken_db      DIR   Kraken2 database path"
    echo "  --bowtie2_index  PATH  Bowtie2 host index prefix"
    echo "  --diamond_db     FILE  DIAMOND database (.dmnd)"
    echo "  --genomad_db     DIR   geNomad database path"
    echo "  --nr_path        FILE  NR protein FASTA path"
    echo ""
    echo "Optional:"
    echo "  --threads        INT   CPU threads (default: 32)"
    echo "  --min_length     INT   Minimum contig length in bp (default: 200)"
    echo "  --skip_existing        Skip steps whose output already exists"
    echo "  --help                 Show this help"
    echo ""
    echo "Example:"
    echo "  bash run_pipeline.sh \\"
    echo "    --fastq_dir     /data/fastq \\"
    echo "    --output_dir    /data/output \\"
    echo "    --kraken_db     /db/kraken \\"
    echo "    --bowtie2_index /db/GR38_bt2 \\"
    echo "    --diamond_db    /db/nr.dmnd \\"
    echo "    --genomad_db    /db/genomad_db \\"
    echo "    --threads 32 --skip_existing"
    exit 1
}

# ── parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --fastq_dir)      FASTQ_DIR="$2";      shift 2 ;;
        --output_dir)     OUTPUT_DIR="$2";     shift 2 ;;
        --kraken_db)      KRAKEN_DB="$2";      shift 2 ;;
        --bowtie2_index)  BOWTIE2_INDEX="$2";  shift 2 ;;
        --diamond_db)     DIAMOND_DB="$2";     shift 2 ;;
        --genomad_db)     GENOMAD_DB="$2";     shift 2 ;;
        --nr_path)        NR_PATH="$2";        shift 2 ;;
        --threads)        THREADS="$2";        shift 2 ;;
        --min_length)     MIN_LENGTH="$2";     shift 2 ;;
        --skip_existing)  SKIP_EXISTING=true;  shift ;;
        --help|-h)        usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# ── validate ──────────────────────────────────────────────────────────────────
MISSING=""
[[ -z "$FASTQ_DIR"     ]] && MISSING+="  --fastq_dir\n"
[[ -z "$OUTPUT_DIR"    ]] && MISSING+="  --output_dir\n"
[[ -z "$KRAKEN_DB"     ]] && MISSING+="  --kraken_db\n"
[[ -z "$BOWTIE2_INDEX" ]] && MISSING+="  --bowtie2_index\n"
[[ -z "$DIAMOND_DB"    ]] && MISSING+="  --diamond_db\n"
[[ -z "$GENOMAD_DB"    ]] && MISSING+="  --genomad_db\n"
[[ -z "$NR_PATH"       ]] && MISSING+="  --nr_path\n"

if [[ -n "$MISSING" ]]; then
    echo "ERROR: Missing required arguments:"
    echo -e "$MISSING"
    usage
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKIP_FLAG=""
[[ "${SKIP_EXISTING}" == "true" ]] && SKIP_FLAG="--skip_existing"

echo "=========================================="
echo "Full pipeline started at $(date)"
echo "Skip existing: ${SKIP_EXISTING}"
echo "=========================================="

bash "${SCRIPT_DIR}/run_phase1.sh" \
    --fastq_dir      "${FASTQ_DIR}" \
    --output_dir     "${OUTPUT_DIR}" \
    --kraken_db      "${KRAKEN_DB}" \
    --bowtie2_index  "${BOWTIE2_INDEX}" \
    --threads        "${THREADS}"

bash "${SCRIPT_DIR}/run_phase2.sh" \
    --fastq_dir   "${FASTQ_DIR}" \
    --output_dir  "${OUTPUT_DIR}" \
    --kraken_db   "${KRAKEN_DB}" \
    --diamond_db  "${DIAMOND_DB}" \
    --genomad_db  "${GENOMAD_DB}" \
    --nr_path     "${NR_PATH}" \
    --threads     "${THREADS}"

bash "${SCRIPT_DIR}/run_phase3.sh" \
    --output_dir  "${OUTPUT_DIR}" \
    --kraken_db   "${KRAKEN_DB}" \
    --diamond_db  "${DIAMOND_DB}" \
    --genomad_db  "${GENOMAD_DB}" \
    --threads     "${THREADS}" \
    --min_length  "${MIN_LENGTH}" \
    ${SKIP_FLAG}

echo ""
echo "=========================================="
echo "Full pipeline finished at $(date)"
echo "=========================================="
