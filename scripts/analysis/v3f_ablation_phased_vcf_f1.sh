#!/bin/bash
# Compare F1 of phased VCFs (baseline vs V5 vs V3F-no-pononly) at 0.6 purity vs SEQC2 truth
# longphase-to phase 不改 FILTER 欄位，預期三版本 F1 完全相同（≈ caller raw F1 = 0.6273）

set -e

OUT=/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation
mkdir -p "$OUT"

TRUTH=/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz
HC_BED=/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed

# 三 phased VCF（0.6 purity）
declare -A VCFS=(
    ["B1_baseline_06"]="/big7_disk/liaoyoyo2001/longphase-to-mod/output/purity_06_simulation/baseline_06/tumor_phased.vcf"
    ["B3_v3f_no_pononly_06"]="/big7_disk/liaoyoyo2001/longphase-to-mod/output/purity_06_simulation/v3f_no_pononly_06/tumor_phased.vcf"
    ["B5_v5_06"]="/big7_disk/liaoyoyo2001/longphase-to-mod/output/purity_06_simulation/v5_06/tumor_phased.vcf"
)

# 0.93 三 phased VCF（對比）
declare -A VCFS_93=(
    ["A1_baseline_093"]="/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_phased.vcf"
    ["A3_v3f_no_pononly_093"]="/big7_disk/liaoyoyo2001/longphase-to-mod/output/v3f_no_pononly/tumor_phased.vcf"
    ["A5_v5_093"]="/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v2b/tumor_phased.vcf"
)

calc_f1() {
    local label=$1
    local vcf=$2

    if [ ! -f "$vcf" ]; then
        echo -e "$label\t-\t-\t-\t-\t-\tMISSING"
        return
    fi

    local TMP=$(mktemp -d)

    # bgzip + index if not gzipped
    local in_vcf="$vcf"
    if [[ "$vcf" != *.gz ]]; then
        bgzip -c "$vcf" > "$TMP/in.vcf.gz" 2>/dev/null
        tabix -f -p vcf "$TMP/in.vcf.gz" 2>/dev/null
        in_vcf="$TMP/in.vcf.gz"
    fi

    bcftools view -f PASS -R "$HC_BED" "$in_vcf" -Oz -o "$TMP/pass.vcf.gz" 2>/dev/null
    tabix -f -p vcf "$TMP/pass.vcf.gz" 2>/dev/null

    bcftools view -R "$HC_BED" "$TRUTH" -Oz -o "$TMP/truth_hc.vcf.gz" 2>/dev/null
    tabix -f -p vcf "$TMP/truth_hc.vcf.gz" 2>/dev/null

    bcftools isec "$TMP/pass.vcf.gz" "$TMP/truth_hc.vcf.gz" -p "$TMP/isec" 2>/dev/null

    local fp=$(grep -v "^#" "$TMP/isec/0000.vcf" | wc -l)
    local fn=$(grep -v "^#" "$TMP/isec/0001.vcf" | wc -l)
    local tp=$(grep -v "^#" "$TMP/isec/0002.vcf" | wc -l)

    local p=$(awk -v tp=$tp -v fp=$fp 'BEGIN{if(tp+fp>0) printf "%.4f", tp/(tp+fp); else print "0"}')
    local r=$(awk -v tp=$tp -v fn=$fn 'BEGIN{if(tp+fn>0) printf "%.4f", tp/(tp+fn); else print "0"}')
    local f1=$(awk -v p=$p -v r=$r 'BEGIN{if(p+r>0) printf "%.4f", 2*p*r/(p+r); else print "0"}')

    echo -e "$label\t$tp\t$fp\t$fn\t$p\t$r\t$f1"

    rm -rf "$TMP"
}

OUT_TSV="$OUT/phased_vcf_f1.tsv"
echo -e "label\tTP\tFP\tFN\tprecision\trecall\tF1" | tee "$OUT_TSV"

# 0.6
for k in B1_baseline_06 B3_v3f_no_pononly_06 B5_v5_06; do
    calc_f1 "$k" "${VCFS[$k]}" | tee -a "$OUT_TSV"
done

# 0.93
for k in A1_baseline_093 A3_v3f_no_pononly_093 A5_v5_093; do
    calc_f1 "$k" "${VCFS_93[$k]}" | tee -a "$OUT_TSV"
done

echo ""
echo "Output: $OUT_TSV"
