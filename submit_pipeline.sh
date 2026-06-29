#!/bin/bash
# Run this ONCE on the LOGIN NODE to copy FASTQs and submit both phases.
# Usage: bash submit_pipeline.sh

RUNDIR=/mnt/hpc_acegid/home/soumareh/260424_VH00635_6_AAF2HNTHV
FASTQ_SRC=${RUNDIR}/fastq                      # original path (symlinks → acegid_storage)
BASE=/mnt/hpc_acegid/nfsscratch/soumareh
FASTQ_DIR=${BASE}/fastq                        # compute-node-accessible copy

mkdir -p "${BASE}/logs" "${FASTQ_DIR}"

# ---- Step 1: copy real FASTQ files to nfsscratch (follows symlinks) ----
echo "Copying FASTQs to nfsscratch (following symlinks)..."
echo "Source: ${FASTQ_SRC}"
echo "Dest:   ${FASTQ_DIR}"
rsync -avL --progress "${FASTQ_SRC}/" "${FASTQ_DIR}/"

# ---- Step 2: count samples ----
N=$(ls "${FASTQ_DIR}"/*_R1*.fastq.gz 2>/dev/null | wc -l)
echo "Found ${N} samples in ${FASTQ_DIR}"

if [[ $N -eq 0 ]]; then
    echo "ERROR: No *_R1*.fastq.gz files found after rsync. Check source path."
    exit 1
fi

echo "First 5 samples:"
ls "${FASTQ_DIR}"/*_R1*.fastq.gz | sort | head -5 \
  | xargs -I{} basename {} | sed 's/_R1.*\.fastq\.gz$//'

# ---- Step 3: submit Phase 1 array (0-indexed) ----
PHASE1_JOB=$(sbatch \
  --array=0-$((N - 1))%8 \
  --parsable \
  phase1_assembly_array.sh)
echo "Phase 1 submitted: job ${PHASE1_JOB}  (${N} tasks, max 8 at once)"

# ---- Step 4: submit Phase 2 — waits for ALL Phase 1 tasks ----
PHASE2_JOB=$(sbatch \
  --dependency=afterok:${PHASE1_JOB} \
  --parsable \
  phase2_diamond.sh)
echo "Phase 2 submitted: job ${PHASE2_JOB}  (starts after Phase 1 completes)"

# ---- Step 5: submit Phase 3 — waits for Phase 2 ----
PHASE3_JOB=$(sbatch \
  --dependency=afterok:${PHASE2_JOB} \
  --parsable \
  phase3_viral_classification.sh)
echo "Phase 3 submitted: job ${PHASE3_JOB}  (starts after Phase 2 completes)"

echo ""
echo "Monitor:"
echo "  squeue -u \$USER"
echo "  tail -f ${BASE}/logs/phase1_${PHASE1_JOB}_0.out"
echo "  tail -f ${BASE}/logs/phase3_${PHASE3_JOB}.out"
