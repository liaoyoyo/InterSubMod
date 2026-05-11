#!/bin/bash
# R1-Global: Phase 2 TO arm — full 40,239 variant global run
# Uses full 144GB Normal BAM via index-based random access (not region-subset).
#
# Paired arm abandoned per user 2026-04-21 (CL-025b concluded abandoned for F1-filter).

set -euo pipefail

BIN="/big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod"
OUT_ROOT="/big7_disk/liaoyoyo2001/InterSubMod/output/hcc1395_normal_pilot_global"
NORMAL_BAM="/big8_disk/liaoyoyo2001/data/bam/HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam"
REF="/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"

# TO arm paths
TO_TUMOR_BAM="/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_tagged.bam"
TO_VCF="$OUT_ROOT/HCC1395_TO_global.vcf.gz"
TO_LOH_BED="/big7_disk/liaoyoyo2001/big7_disk_output/bip8_output_archive/research_rounds/20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_phased_LOH.bed"

THREADS=16
WINDOW=5000

mkdir -p "$OUT_ROOT/TO" "$OUT_ROOT/logs"

echo "[+] Launching TO arm global pipeline at $(date)"
echo "    Variants: 40,239"
echo "    Normal BAM: $NORMAL_BAM (full 144GB, random access)"
echo "    Tumor BAM: $TO_TUMOR_BAM"
echo "    LOH bed: $TO_LOH_BED"
echo "    Output: $OUT_ROOT/TO/"

"$BIN" \
    -t "$TO_TUMOR_BAM" \
    -n "$NORMAL_BAM" \
    --loh-bed "$TO_LOH_BED" \
    -r "$REF" \
    -v "$TO_VCF" \
    -o "$OUT_ROOT/TO" \
    -w $WINDOW \
    -j $THREADS \
    --distance-metric BERNOULLI \
    --log-level info \
    > "$OUT_ROOT/logs/to_global.log" 2>&1

EXIT_CODE=$?
echo "[+] TO global done (exit=$EXIT_CODE) at $(date)"

if [ $EXIT_CODE -eq 0 ]; then
    echo "--- Output summary ---"
    ls -la "$OUT_ROOT/TO/" | head -20
    echo "--- Log tail ---"
    tail -20 "$OUT_ROOT/logs/to_global.log"
else
    echo "[-] FAILED"
    tail -40 "$OUT_ROOT/logs/to_global.log"
    exit 1
fi
