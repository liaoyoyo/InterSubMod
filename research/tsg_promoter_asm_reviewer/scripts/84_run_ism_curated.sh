#!/usr/bin/env bash
# 84 - Re-run ISM on the curated locus set (HCC1395 TP + FP) WITH distance-matrix
# output enabled (existence scan used --no-output-distance-matrix). Produces the
# per-region distance/NHD/matrix.csv + methylation/methylation.csv + reads/reads.tsv
# the redesigned display renders both heatmaps from.
set -eu
export TMPDIR=/big7_disk/liaoyoyo2001/tmp; mkdir -p "$TMPDIR"
ISM=/big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod
REF=/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta
LPS=/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s
TBAM="$LPS/HCC1395_tagged.bam"
DV=/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2/cn_confound/cross_sample/display_v2
OUTROOT=/big7_disk/liaoyoyo2001/ism_display_scan
THREADS=${THREADS:-16}
mkdir -p "$OUTROOT"

run_cls() {
  local cls=$1
  local src="$LPS/filtered_snv_${cls}.vcf.gz"
  local sub="$DV/curated_${cls}.vcf.gz"
  local out="$OUTROOT/HCC1395_${cls}"
  echo "[84] $cls: subset VCF from $(basename "$src")"
  bcftools view -T "$DV/${cls}_pos.txt" "$src" -Oz -o "$sub"
  tabix -f -p vcf "$sub"
  local nrec; nrec=$(zcat "$sub" | grep -vc '^#')
  echo "[84] $cls: $nrec records -> running ISM (distance matrix ON)"
  rm -rf "$out"; mkdir -p "$out"
  local t0=$SECONDS
  "$ISM" --tumor-bam "$TBAM" --reference "$REF" --vcf "$sub" \
         --output-dir "$out" --threads "$THREADS" --window-size 1000 \
         --distance-metric NHD --log-level error > "$out/run.log" 2>&1
  local rc=$?
  echo "[84] $cls: rc=$rc elapsed=$((SECONDS-t0))s"
  # prune bulk we don't render (strand matrices, forward/reverse meth, debug)
  find "$out" -type f \( -name 'matrix_forward.csv' -o -name 'matrix_reverse.csv' \
       -o -name 'methylation_forward.csv' -o -name 'methylation_reverse.csv' \
       -o -name 'stats_forward.txt' -o -name 'stats_reverse.txt' \
       -o -name 'tree_forward.nwk' -o -name 'tree_reverse.nwk' \) -delete 2>/dev/null || true
  rm -rf "$out/debug" 2>/dev/null || true
  local n; n=$( [ -f "$out/significance_summary.csv" ] && echo $(( $(wc -l < "$out/significance_summary.csv") - 1 )) || echo NA )
  echo "[84] $cls: summary rows=$n  dir size=$(du -sh "$out" 2>/dev/null | cut -f1)"
}

run_cls tp
run_cls fp
echo "[84] ===== CURATED ISM DISPLAY SCAN DONE ====="
for cls in tp fp; do
  echo "  HCC1395_${cls}: $(find "$OUTROOT/HCC1395_${cls}" -name 'matrix.csv' | wc -l) regions with distance matrix"
done
