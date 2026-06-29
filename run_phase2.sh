#!/bin/bash
# Phase 2: collect contigs, run geNomad + DIAMOND + alignment
#
# Usage:
#   bash run_phase2.sh --fastq_dir <path> --output_dir <path> \
#                      --kraken_db <path> --diamond_db <path> \
#                      --genomad_db <path> --nr_path <path> [--threads 32]

set -eo pipefail

# ── defaults ─────────────────────────────────────────────────────────────────
THREADS=32

usage() {
    echo ""
    echo "Usage: bash run_phase2.sh [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  --fastq_dir    DIR   Directory with *_R1*.fastq.gz files (for sample discovery)"
    echo "  --output_dir   DIR   Output directory (must contain per-sample subdirs from phase 1)"
    echo "  --kraken_db    DIR   Kraken2 database path"
    echo "  --diamond_db   FILE  DIAMOND database (.dmnd)"
    echo "  --genomad_db   DIR   geNomad database path"
    echo "  --nr_path      FILE  NR protein FASTA path"
    echo ""
    echo "Optional:"
    echo "  --threads      INT   CPU threads (default: 32)"
    echo "  --help               Show this help"
    echo ""
    echo "Example:"
    echo "  bash run_phase2.sh \\"
    echo "    --fastq_dir   /data/fastq \\"
    echo "    --output_dir  /data/output \\"
    echo "    --kraken_db   /db/kraken \\"
    echo "    --diamond_db  /db/nr_genomad.dmnd \\"
    echo "    --genomad_db  /db/genomad_db \\"
    echo "    --nr_path     /db/nr \\"
    echo "    --threads 32"
    exit 1
}

# ── parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --fastq_dir)   FASTQ_DIR="$2";   shift 2 ;;
        --output_dir)  OUTPUT_DIR="$2";  shift 2 ;;
        --kraken_db)   KRAKEN_DB="$2";   shift 2 ;;
        --diamond_db)  DIAMOND_DB="$2";  shift 2 ;;
        --genomad_db)  GENOMAD_DB="$2";  shift 2 ;;
        --nr_path)     NR_PATH="$2";     shift 2 ;;
        --threads)     THREADS="$2";     shift 2 ;;
        --help|-h)     usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# ── validate ──────────────────────────────────────────────────────────────────
MISSING=""
[[ -z "$FASTQ_DIR"  ]] && MISSING+="  --fastq_dir\n"
[[ -z "$OUTPUT_DIR" ]] && MISSING+="  --output_dir\n"
[[ -z "$KRAKEN_DB"  ]] && MISSING+="  --kraken_db\n"
[[ -z "$DIAMOND_DB" ]] && MISSING+="  --diamond_db\n"
[[ -z "$GENOMAD_DB" ]] && MISSING+="  --genomad_db\n"
[[ -z "$NR_PATH"    ]] && MISSING+="  --nr_path\n"

if [[ -n "$MISSING" ]]; then
    echo "ERROR: Missing required arguments:"
    echo -e "$MISSING"
    usage
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH}"
mkdir -p "${OUTPUT_DIR}/logs"

echo "=========================================="
echo "Phase 2 started at $(date)"
echo "Output dir : ${OUTPUT_DIR}"
echo "=========================================="

# ── collect contigs and kraken reports into flat output dir ───────────────────
echo "Collecting contigs from per-sample directories..."
for SAMPLE_OUT in "${OUTPUT_DIR}"/*/; do
    SAMPLE=$(basename "${SAMPLE_OUT}")
    SRC="${SAMPLE_OUT}/${SAMPLE}_contigs.fasta"
    DST="${OUTPUT_DIR}/${SAMPLE}_contigs.fasta"
    [[ -f "$SRC" && ! -f "$DST" ]] && cp "$SRC" "$DST" && echo "  Collected: ${SAMPLE}_contigs.fasta"

    SRC="${SAMPLE_OUT}/${SAMPLE}_kraken_report.txt"
    DST="${OUTPUT_DIR}/${SAMPLE}_kraken_report.txt"
    [[ -f "$SRC" && ! -f "$DST" ]] && cp "$SRC" "$DST"
done

N=$(ls "${OUTPUT_DIR}"/*_contigs.fasta 2>/dev/null | wc -l)
echo "Total contig files: ${N}"

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
    2>&1 | tee "${OUTPUT_DIR}/logs/phase2.log"

echo ""
echo "=========================================="
echo "Phase 2 complete at $(date)"
echo "=========================================="
