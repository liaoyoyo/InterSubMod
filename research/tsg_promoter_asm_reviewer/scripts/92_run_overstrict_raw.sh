#!/usr/bin/env bash
# 92 - Re-run ISM on ALL over-strict loci (clean HP-aligned PERMANOVA structure but
# gated CramersV<0.1) to obtain RAW CramersV (clustering/significance.json), confirming
# the reliability-gated cause genome-wide (beyond the 425 curated subset). Distance
# matrix ON only to emit per-region significance.json; pruned to keep just that.
set -eu
export TMPDIR=/big7_disk/liaoyoyo2001/tmp; mkdir -p "$TMPDIR"
ISM=/big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod
REF=/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta
LPS=/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s
TBAM="$LPS/HCC1395_tagged.bam"
DV=/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2/cn_confound/cross_sample/display_v2
OUTROOT=/big7_disk/liaoyoyo2001/ism_overstrict_raw
THREADS=${THREADS:-16}
mkdir -p "$OUTROOT"

run_cls() {
  local cls=$1
  local sub="$DV/overstrict_${cls}.vcf.gz"
  local out="$OUTROOT/HCC1395_${cls}"
  bcftools view -T "$DV/overstrict_${cls}_pos.txt" "$LPS/filtered_snv_${cls}.vcf.gz" -Oz -o "$sub"
  tabix -f -p vcf "$sub"
  local nrec; nrec=$(zcat "$sub" | grep -vc '^#')
  echo "[92] $cls: $nrec records -> ISM (raw CramersV via significance.json)"
  rm -rf "$out"; mkdir -p "$out"
  local t0=$SECONDS
  "$ISM" --tumor-bam "$TBAM" --reference "$REF" --vcf "$sub" \
         --output-dir "$out" --threads "$THREADS" --window-size 1000 \
         --distance-metric NHD --log-level error > "$out/run.log" 2>&1
  echo "[92] $cls: rc=$? elapsed=$((SECONDS-t0))s"
  # keep ONLY clustering/significance.json per region; drop bulky distance/methylation/reads
  find "$out" -type d \( -name distance -o -name methylation -o -name reads \) -exec rm -rf {} + 2>/dev/null || true
  find "$out" -type f \( -name 'leaf_order.txt' -o -name 'linkage_matrix.csv' -o -name 'tree*.nwk' \) -delete 2>/dev/null || true
  rm -rf "$out/debug" 2>/dev/null || true
  echo "[92] $cls: kept $(find "$out" -name significance.json | wc -l) significance.json; size=$(du -sh "$out" | cut -f1)"
}

run_cls tp
run_cls fp
echo "[92] ===== OVER-STRICT RAW RE-RUN DONE ====="
