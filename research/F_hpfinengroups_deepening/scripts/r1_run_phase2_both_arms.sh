#!/bin/bash
# R1 HCC1395 Normal BAM pilot — Phase 2 execution (TO + paired in parallel)
#
# TO arm:     tumor_tagged.bam + HCC1395_normal_subset.bam + TO LOH bed + TO subset VCF
# Paired arm: HCC1395_Tmode_tagged.bam + HCC1395_normal_subset.bam + paired subset VCF
#             (paired mode uses ClairS Paired internal LOH, no --loh-bed flag)
#
# Both outputs go to output/hcc1395_normal_pilot/{TO,paired}/

set -euo pipefail

BIN="/big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod"
OUT_ROOT="/big7_disk/liaoyoyo2001/InterSubMod/output/hcc1395_normal_pilot"
NORMAL_BAM="$OUT_ROOT/HCC1395_normal_subset.bam"
REF="/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"

# TO arm paths
TO_TUMOR_BAM="/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_tagged.bam"
TO_VCF="$OUT_ROOT/HCC1395_TO_filter_subset.vcf.gz"
TO_LOH_BED="/big7_disk/liaoyoyo2001/big7_disk_output/bip8_output_archive/research_rounds/20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_phased_LOH.bed"

# Paired arm paths
PAIRED_TUMOR_BAM="/big8_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam"
PAIRED_VCF="$OUT_ROOT/HCC1395_paired_filter_subset.vcf.gz"

THREADS=16
WINDOW=5000

mkdir -p "$OUT_ROOT/TO" "$OUT_ROOT/paired"
mkdir -p "$OUT_ROOT/logs"

launch_to() {
    echo "[+] Launching TO arm pipeline..."
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
        > "$OUT_ROOT/logs/to_arm.log" 2>&1
    echo "[+] TO arm done (exit=$?)"
}

launch_paired() {
    echo "[+] Launching paired arm pipeline..."
    "$BIN" \
        -t "$PAIRED_TUMOR_BAM" \
        -n "$NORMAL_BAM" \
        -r "$REF" \
        -v "$PAIRED_VCF" \
        -o "$OUT_ROOT/paired" \
        -w $WINDOW \
        -j $THREADS \
        --distance-metric BERNOULLI \
        --log-level info \
        > "$OUT_ROOT/logs/paired_arm.log" 2>&1
    echo "[+] Paired arm done (exit=$?)"
}

# Run both arms in parallel
launch_to &
TO_PID=$!
launch_paired &
PAIRED_PID=$!

echo "[+] TO PID=$TO_PID, Paired PID=$PAIRED_PID"
echo "[+] Waiting..."

wait $TO_PID
TO_EXIT=$?
wait $PAIRED_PID
PAIRED_EXIT=$?

echo "[+] TO exit=$TO_EXIT, Paired exit=$PAIRED_EXIT"

if [ $TO_EXIT -eq 0 ] && [ $PAIRED_EXIT -eq 0 ]; then
    echo "[+] Both arms SUCCESS"
    ls -la "$OUT_ROOT/TO/" | head -20
    echo "---"
    ls -la "$OUT_ROOT/paired/" | head -20
    exit 0
else
    echo "[-] At least one arm FAILED"
    echo "--- TO arm tail ---"
    tail -30 "$OUT_ROOT/logs/to_arm.log"
    echo "--- Paired arm tail ---"
    tail -30 "$OUT_ROOT/logs/paired_arm.log"
    exit 1
fi
