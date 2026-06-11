#!/usr/bin/env bash
# 71 - Complete ISM existence scan: 6 samples × {TP,FP,FN} via unmodified ISM.
# Native Stage-2 HP-fine + per-CpG Fisher + signed Δβ + Significant gate.
# No normal (existence = tumor-only HP-axis). significance_summary.csv per (sample,class).
set -u
export TMPDIR=/big7_disk/liaoyoyo2001/tmp; mkdir -p "$TMPDIR"
ISM=/big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod
REF=/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta
CANON=/big7_disk/liaoyoyo2001/big7_disk_output/canonical
OUTROOT=/big7_disk/liaoyoyo2001/ism_existence_scan
THREADS=${THREADS:-8}
mkdir -p "$OUTROOT"

run_one() {
  local s=$1 cls=$2 rundir=$3
  local lps="$CANON/$s/paired_full/$rundir/longphase_s"
  local tbam="$lps/${s}_tagged.bam"
  local vcf="$lps/filtered_snv_${cls}.vcf.gz"
  local out="$OUTROOT/${s}_${cls}"
  [ -f "$tbam" ] || { echo "[71] MISS tumor $tbam"; return; }
  [ -f "$vcf" ] || { echo "[71] MISS vcf $vcf"; return; }
  mkdir -p "$out"
  local t0=$SECONDS
  "$ISM" --tumor-bam "$tbam" --reference "$REF" --vcf "$vcf" \
         --output-dir "$out" --threads "$THREADS" --window-size 1000 \
         --no-output-distance-matrix --log-level error > "$out/run.log" 2>&1
  local rc=$?
  local n=$( [ -f "$out/significance_summary.csv" ] && echo $(( $(wc -l < "$out/significance_summary.csv") - 1 )) || echo NA )
  echo "[71] $s $cls rc=$rc rows=$n elapsed=$((SECONDS-t0))s"
  # prune bulky per-region/debug output, keep the summary CSVs
  rm -rf "$out/debug" "$out/filtered_snv_${cls}" 2>/dev/null
}

# per-sample stream (tp,fp,fn sequential); 6 samples in parallel
declare -A RUN=( [HCC1395]=20260314_HCC1395_paired_full_full_complete_matrix )
for s in HCC1937 HCC1954 H1437 H2009 COLO829; do RUN[$s]=20260315_${s}_paired_full_full_complete_matrix; done

for s in HCC1395 HCC1937 HCC1954 H1437 H2009 COLO829; do
  ( for cls in tp fp fn; do run_one "$s" "$cls" "${RUN[$s]}"; done ) &
done
wait
echo "[71] ===== ALL ISM EXISTENCE SCANS DONE ====="
for s in HCC1395 HCC1937 HCC1954 H1437 H2009 COLO829; do
  for cls in tp fp fn; do
    f="$OUTROOT/${s}_${cls}/significance_summary.csv"
    [ -f "$f" ] && echo "  ${s}_${cls}: $(( $(wc -l < "$f") - 1 )) rows"
  done
done
