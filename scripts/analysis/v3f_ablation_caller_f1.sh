#!/bin/bash
# V3F Ablation Caller F1 vs SEQC2 truth
# 計算 ClairS-TO PASS variants vs SEQC2 high-confidence sSNV truth
# 對 0.93 純樣本 + 0.6 t30_n20 兩種 purity

set -e

OUT=/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation
mkdir -p "$OUT"

TRUTH=/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz
HC_BED=/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed

VCF_093=/big7_disk/liaoyoyo2001/data/HCC1395_5kHz/snv.vcf.gz
VCF_06=/big8_disk/data/HCC1395/ONT/subsample/t30_n20/ClairS_TO_v0_3_0/snv.vcf.gz

calc_f1() {
    local label=$1
    local clairs_vcf=$2

    local TMP=$(mktemp -d)

    # 1. PASS filter + restrict to HC regions
    bcftools view -f PASS -R "$HC_BED" "$clairs_vcf" -Oz -o "$TMP/clairs_pass.vcf.gz" 2>/dev/null
    tabix -f -p vcf "$TMP/clairs_pass.vcf.gz"

    # 2. Truth restrict to HC regions (already in HC, 但確認)
    bcftools view -R "$HC_BED" "$TRUTH" -Oz -o "$TMP/truth_hc.vcf.gz" 2>/dev/null
    tabix -f -p vcf "$TMP/truth_hc.vcf.gz"

    # 3. Intersect
    bcftools isec "$TMP/clairs_pass.vcf.gz" "$TMP/truth_hc.vcf.gz" -p "$TMP/isec" 2>/dev/null

    # 4. Count
    local fp=$(grep -v "^#" "$TMP/isec/0000.vcf" | wc -l)
    local fn=$(grep -v "^#" "$TMP/isec/0001.vcf" | wc -l)
    local tp=$(grep -v "^#" "$TMP/isec/0002.vcf" | wc -l)

    # 5. F1
    local precision=$(awk -v tp=$tp -v fp=$fp 'BEGIN{printf "%.4f", tp/(tp+fp)}')
    local recall=$(awk -v tp=$tp -v fn=$fn 'BEGIN{printf "%.4f", tp/(tp+fn)}')
    local f1=$(awk -v p=$precision -v r=$recall 'BEGIN{printf "%.4f", 2*p*r/(p+r)}')

    echo -e "$label\t$tp\t$fp\t$fn\t$precision\t$recall\t$f1"

    # Clean
    rm -rf "$TMP"
}

# Header
echo -e "label\tTP\tFP\tFN\tprecision\trecall\tF1" | tee "$OUT/caller_f1.tsv"

# 0.93 純樣本 ClairS-TO PASS F1
calc_f1 "ClairS-TO_0.93_raw" "$VCF_093" | tee -a "$OUT/caller_f1.tsv"

# 0.6 t30_n20 ClairS-TO PASS F1
calc_f1 "ClairS-TO_0.6_raw" "$VCF_06" | tee -a "$OUT/caller_f1.tsv"

echo ""
echo "Output: $OUT/caller_f1.tsv"
