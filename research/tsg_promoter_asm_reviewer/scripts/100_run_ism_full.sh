#!/usr/bin/env bash
# 100 - Re-run ISM on the FULL HCC1395 TP + FP somatic SNVs WITH distance-matrix output,
# so every one of the 30,350 ISM-analysed loci gets per-region methylation.csv +
# distance/NHD/matrix.csv + reads.tsv to render BOTH heatmaps. Classification stays on
# the frozen existence-scan summary (verified identical); this run only supplies figures.
# Prune to render essentials to bound disk (~2.4 GB).
set -eu
export TMPDIR=/big7_disk/liaoyoyo2001/tmp; mkdir -p "$TMPDIR"
ISM=/big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod
REF=/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta
LPS=/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s
TBAM="$LPS/HCC1395_tagged.bam"
OUTROOT=/big7_disk/liaoyoyo2001/ism_full_scan
THREADS=${THREADS:-24}
mkdir -p "$OUTROOT"

run_cls() {
  local cls=$1
  local vcf="$LPS/filtered_snv_${cls}.vcf.gz"
  local out="$OUTROOT/HCC1395_${cls}"
  local nrec; nrec=$(zcat "$vcf" | grep -vc '^#')
  echo "[100] $cls: $nrec records -> ISM (full, distance matrix ON)"
  rm -rf "$out"; mkdir -p "$out"
  local t0=$SECONDS
  "$ISM" --tumor-bam "$TBAM" --reference "$REF" --vcf "$vcf" \
         --output-dir "$out" --threads "$THREADS" --window-size 1000 \
         --distance-metric NHD --log-level error > "$out/run.log" 2>&1
  echo "[100] $cls: rc=$? elapsed=$((SECONDS-t0))s"
  # keep ONLY what the renderer needs: methylation.csv, distance/NHD/matrix.csv, reads.tsv,
  # clustering/{leaf_order.txt,significance.json}. Drop strand matrices / fwd-rev meth / trees / debug.
  find "$out" -type f \( -name 'matrix_forward.csv' -o -name 'matrix_reverse.csv' \
       -o -name 'methylation_forward.csv' -o -name 'methylation_reverse.csv' \
       -o -name 'stats_forward.txt' -o -name 'stats_reverse.txt' \
       -o -name 'tree*.nwk' -o -name 'linkage_matrix.csv' \) -delete 2>/dev/null || true
  rm -rf "$out/debug" 2>/dev/null || true
  local n; n=$(find "$out" -name 'matrix.csv' | wc -l)
  echo "[100] $cls: regions with distance matrix=$n  size=$(du -sh "$out" 2>/dev/null | cut -f1)  disk_free=$(df -h /big7_disk | tail -1 | awk '{print $4}')"
}

run_cls tp
run_cls fp
echo "[100] ===== FULL ISM SCAN DONE ====="
for cls in tp fp; do echo "  HCC1395_${cls}: $(find "$OUTROOT/HCC1395_${cls}" -name 'matrix.csv' | wc -l) regions"; done
