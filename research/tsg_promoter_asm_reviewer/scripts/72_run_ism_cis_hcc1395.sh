#!/usr/bin/env bash
# 72 - HCC1395 normal-anchored cis-ASM scan: TP/FP/FN with tumor + FULL normal.
# ISM native HP_Residual_Delta = tumor_HP_delta - normal_HP_delta = cis signal
# (somatic-allele methylation change vs the matched-normal baseline). Plus SampleASM_*.
set -u
export TMPDIR=/big7_disk/liaoyoyo2001/tmp; mkdir -p "$TMPDIR"
ISM=/big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod
REF=/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta
LPS=/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s
TBAM="$LPS/HCC1395_tagged.bam"
NBAM=/big8_disk/liaoyoyo2001/data/bam/HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam
OUTROOT=/big7_disk/liaoyoyo2001/ism_cis_scan
THREADS=${THREADS:-32}
mkdir -p "$OUTROOT"

for cls in tp fp fn; do
  vcf="$LPS/filtered_snv_${cls}.vcf.gz"
  out="$OUTROOT/HCC1395_${cls}_cis"
  [ -f "$vcf" ] || { echo "[72] MISS $vcf"; continue; }
  mkdir -p "$out"
  t0=$SECONDS
  "$ISM" --tumor-bam "$TBAM" --normal-bam "$NBAM" --reference "$REF" --vcf "$vcf" \
         --output-dir "$out" --threads "$THREADS" --window-size 1000 \
         --no-output-distance-matrix --log-level error > "$out/run.log" 2>&1
  rc=$?
  n=$( [ -f "$out/significance_summary.csv" ] && echo $(( $(wc -l < "$out/significance_summary.csv") - 1 )) || echo NA )
  echo "[72] HCC1395 $cls cis rc=$rc rows=$n elapsed=$((SECONDS-t0))s"
  rm -rf "$out/debug" "$out/filtered_snv_${cls}" 2>/dev/null
done
echo "[72] ===== HCC1395 CIS SCAN DONE ====="
