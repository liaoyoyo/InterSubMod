#!/usr/bin/env bash
# Launch BAM-QC extraction across all (sample, mode) pairs in parallel.
#
# Strategy:
#   - 4 concurrent BAMs, each using 4 pysam I/O threads (16 threads total)
#   - Outputs are TSV per combo under research/feature_layered_observation/data/bam_readqc/
#   - Idempotent: skips combos whose output already exists (non-empty)
#   - Progress logged to research/feature_layered_observation/data/bam_readqc_progress.log
#
# Usage:
#     nohup bash research/feature_layered_observation/scripts/run_bam_readqc_all.sh \
#         > research/feature_layered_observation/data/bam_readqc_nohup.out 2>&1 &

set -u

ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
SCRIPT="$ROOT/research/feature_layered_observation/scripts/bam_readqc_per_region.py"
LOG="$ROOT/research/feature_layered_observation/data/bam_readqc_progress.log"
OUT_DIR="$ROOT/research/feature_layered_observation/data/bam_readqc"
mkdir -p "$OUT_DIR"

# (sample, mode) pairs matching master_extended coverage.
# COLO829 to_pileup is included — BAM exists but no master_extended rows;
# the script will produce an empty output (0 regions) which is harmless.
PAIRS=(
    "HCC1395 paired_full"
    "HCC1395_DORADO paired_full"
    "H1437 paired_full"
    "H2009 paired_full"
    "HCC1937 paired_full"
    "HCC1954 paired_full"
    "COLO829 paired_full"
    "HCC1395 to_pileup"
    "HCC1395_DORADO to_pileup"
    "H1437 to_pileup"
    "H2009 to_pileup"
    "HCC1937 to_pileup"
    "HCC1954 to_pileup"
    "COLO829 to_pileup"
)

THREADS_PER_BAM=4
MAX_PARALLEL=4

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [launch] pairs=${#PAIRS[@]} max_parallel=$MAX_PARALLEL threads/bam=$THREADS_PER_BAM" | tee -a "$LOG"

# xargs-based parallel pool
printf '%s\n' "${PAIRS[@]}" | \
xargs -d '\n' -P "$MAX_PARALLEL" -I {} bash -c '
    pair="{}"
    sample=$(echo "$pair" | awk "{print \$1}")
    mode=$(echo "$pair" | awk "{print \$2}")
    out="'"$OUT_DIR"'/bam_readqc_${sample}_${mode}.tsv"
    if [ -s "$out" ]; then
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] [skip] $sample/$mode already has $out" >> "'"$LOG"'"
        exit 0
    fi
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] [dispatch] $sample/$mode" >> "'"$LOG"'"
    python3 "'"$SCRIPT"'" \
        --sample "$sample" --mode "$mode" \
        --threads '"$THREADS_PER_BAM"' \
        --log "'"$LOG"'" \
        >> "'"$LOG"'" 2>&1
    rc=$?
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] [return] $sample/$mode rc=$rc" >> "'"$LOG"'"
'

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [all-done] driver exiting" | tee -a "$LOG"
