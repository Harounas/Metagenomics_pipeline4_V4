#!/bin/bash
# Master script: run all 3 phases in sequence
# Usage: bash run_pipeline.sh
#
# To run a single phase:
#   bash run_phase1.sh
#   bash run_phase2.sh
#   bash run_phase3.sh

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=========================================="
echo "Full pipeline started at $(date)"
echo "=========================================="

bash "${SCRIPT_DIR}/run_phase1.sh"
bash "${SCRIPT_DIR}/run_phase2.sh"
bash "${SCRIPT_DIR}/run_phase3.sh"

echo ""
echo "=========================================="
echo "Full pipeline finished at $(date)"
echo "=========================================="
