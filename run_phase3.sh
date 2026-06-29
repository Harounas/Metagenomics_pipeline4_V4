#!/bin/bash
# Phase 3: Kraken2 + DIAMOND on merged contigs -> CD-HIT cluster -> final viral TSV
#
# Usage:
#   bash run_phase3.sh --output_dir <path> --kraken_db <path> \
#                      --diamond_db <path> --nr_path <path> [--threads 32] [--min_length 200]

set -eo pipefail

# ── defaults ─────────────────────────────────────────────────────────────────
THREADS=32
MIN_LENGTH=200

usage() {
    echo ""
    echo "Usage: bash run_phase3.sh [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  --output_dir   DIR   Directory with *_contigs.fasta files (phase 1 output)"
    echo "  --kraken_db    DIR   Kraken2 database path"
    echo "  --diamond_db   FILE  DIAMOND database (.dmnd)"
    echo "  --nr_path      FILE  NR protein FASTA path"
    echo ""
    echo "Optional:"
    echo "  --threads      INT   CPU threads (default: 32)"
    echo "  --min_length   INT   Minimum contig length in bp (default: 200)"
    echo "  --help               Show this help"
    echo ""
    echo "Example:"
    echo "  bash run_phase3.sh \\"
    echo "    --output_dir  /data/output \\"
    echo "    --kraken_db   /db/kraken \\"
    echo "    --diamond_db  /db/nr_genomad.dmnd \\"
    echo "    --nr_path     /db/nr \\"
    echo "    --threads 32"
    exit 1
}

# ── parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --output_dir)  OUTPUT_DIR="$2";  shift 2 ;;
        --kraken_db)   KRAKEN_DB="$2";   shift 2 ;;
        --diamond_db)  DIAMOND_DB="$2";  shift 2 ;;
        --nr_path)     NR_PATH="$2";     shift 2 ;;
        --threads)     THREADS="$2";     shift 2 ;;
        --min_length)  MIN_LENGTH="$2";  shift 2 ;;
        --help|-h)     usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# ── validate ──────────────────────────────────────────────────────────────────
MISSING=""
[[ -z "$OUTPUT_DIR" ]] && MISSING+="  --output_dir\n"
[[ -z "$KRAKEN_DB"  ]] && MISSING+="  --kraken_db\n"
[[ -z "$DIAMOND_DB" ]] && MISSING+="  --diamond_db\n"
[[ -z "$NR_PATH"    ]] && MISSING+="  --nr_path\n"

if [[ -n "$MISSING" ]]; then
    echo "ERROR: Missing required arguments:"
    echo -e "$MISSING"
    usage
fi

[[ ! -d "$OUTPUT_DIR" ]] && { echo "ERROR: --output_dir not found: $OUTPUT_DIR"; exit 1; }
[[ ! -d "$KRAKEN_DB"  ]] && { echo "ERROR: --kraken_db not found: $KRAKEN_DB";  exit 1; }
[[ ! -f "$DIAMOND_DB" ]] && { echo "ERROR: --diamond_db not found: $DIAMOND_DB"; exit 1; }
[[ ! -f "$NR_PATH"    ]] && { echo "ERROR: --nr_path not found: $NR_PATH";       exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH}"
mkdir -p "${OUTPUT_DIR}/logs"

echo "=========================================="
echo "Phase 3 started at $(date)"
echo "Output dir  : ${OUTPUT_DIR}"
echo "Kraken DB   : ${KRAKEN_DB}"
echo "DIAMOND DB  : ${DIAMOND_DB}"
echo "NR path     : ${NR_PATH}"
echo "Threads     : ${THREADS}"
echo "Min length  : ${MIN_LENGTH} bp"
echo "=========================================="

python "${SCRIPT_DIR}/Metagenomics_pipeline4_V2/viral_classification_workflow.py" \
    --output_dir     "${OUTPUT_DIR}" \
    --kraken_db      "${KRAKEN_DB}" \
    --diamond_db     "${DIAMOND_DB}" \
    --nr_path        "${NR_PATH}" \
    --threads        "${THREADS}" \
    --min_length     "${MIN_LENGTH}" \
    --skip_existing \
    2>&1 | tee "${OUTPUT_DIR}/logs/phase3.log"

FINAL_TSV="${OUTPUT_DIR}/filtered_clusters_assigned_rep_virus.tsv"
if [[ -f "${FINAL_TSV}" ]]; then
    N=$(tail -n +2 "${FINAL_TSV}" | wc -l)
    echo ""
    echo "=========================================="
    echo "Phase 3 complete at $(date)"
    echo "Final TSV : ${FINAL_TSV}"
    echo "Viral rows: ${N}"
    echo "=========================================="
else
    echo "ERROR: Final TSV not found at ${FINAL_TSV}" >&2
    exit 1
fi
