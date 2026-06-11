#!/bin/bash
# Phase C with significance: HCC1395 全基因組 V3F vs V5 vs V6 three-way
# 12 runs (3 BAM × 2 flag × 2 label)
# 2026-05-18 amendment of run_phaseC_genome_three_way.sh:
#   - Removed --no-distance-matrix flag (line 48+61) to enable significance computation
#   - Changed OUT_BASE to phaseC_genome_three_way_with_significance/ (avoid overwriting original)
#   - Plan: InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/00_PLAN.md §Step -1
set -euo pipefail

export TMPDIR=/big7_disk/liaoyoyo2001/tmp
mkdir -p $TMPDIR

REPO=/big7_disk/liaoyoyo2001/InterSubMod
ISM=$REPO/build/bin/inter_sub_mod
REF=/big7_disk/liaoyoyo2001/data/ref/GRCh38_no_alt_analysis_set.fasta
LOH_BED=/big7_disk/liaoyoyo2001/longphase-to-mod/output/threshold_compare/v5_flag/tumor_phased_LOH.bed
VCF_TP=/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz
VCF_FP=/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_fp.vcf.gz
THREADS=24

declare -A BAMS
BAMS[V3F]=/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v3_fixed/tumor_tagged.bam
BAMS[V5]=/big7_disk/liaoyoyo2001/longphase-to-mod/output/threshold_compare/v5_flag/tumor_tagged.bam
BAMS[V6]=/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_germline_absent_revert/tumor_tagged.bam

# 2026-05-18: new OUT_BASE to avoid overwriting original phaseC results
OUT_BASE=$REPO/research/paired_priority_bug_audit/phaseC_genome_three_way_with_significance
mkdir -p "$OUT_BASE"

echo "[Phase C with significance] $(date)"
echo "  TP VCF: $VCF_TP ($(zcat $VCF_TP | grep -v '^#' | wc -l) SNVs)"
echo "  FP VCF: $VCF_FP ($(zcat $VCF_FP | grep -v '^#' | wc -l) SNVs)"

run_ism() {
    local bam_tag=$1
    local flag_label=$2
    local flag_arg=$3
    local vcf_path=$4
    local label=$5
    local out_dir=$OUT_BASE/${bam_tag}_${flag_label}_${label}
    mkdir -p "$out_dir"
    echo "  [run] BAM=$bam_tag flag=$flag_label label=$label → $out_dir"
    local t0=$(date +%s)
    # 2026-05-18: removed --no-distance-matrix to enable significance
    if [ -n "$flag_arg" ]; then
        $ISM \
            --tumor-bam "${BAMS[$bam_tag]}" \
            --reference "$REF" \
            --vcf "$vcf_path" \
            --loh-bed "$LOH_BED" \
            --output-dir "$out_dir" \
            --window-size 5000 \
            --threads $THREADS \
            --log-level warn \
            $flag_arg \
            > "$out_dir/run.log" 2>&1
    else
        $ISM \
            --tumor-bam "${BAMS[$bam_tag]}" \
            --reference "$REF" \
            --vcf "$vcf_path" \
            --loh-bed "$LOH_BED" \
            --output-dir "$out_dir" \
            --window-size 5000 \
            --threads $THREADS \
            --log-level warn \
            > "$out_dir/run.log" 2>&1
    fi
    local t1=$(date +%s)
    echo "    elapsed: $((t1-t0))s"
}

# Allow running specific subset via $1 (e.g. "V6_off_tp" for pilot)
if [ $# -gt 0 ]; then
    case "$1" in
        V*_*_*)
            IFS='_' read -ra parts <<< "$1"
            bam_tag="${parts[0]}"
            flag_label="${parts[1]}"
            label="${parts[2]}"
            vcf_var="VCF_${label^^}"
            vcf_path="${!vcf_var}"
            flag_arg=""
            if [ "$flag_label" = "on" ]; then flag_arg="--germline-hp-only"; fi
            echo "=== Single run: $1 ==="
            run_ism "$bam_tag" "$flag_label" "$flag_arg" "$vcf_path" "$label"
            echo
            echo "[Phase C single $1] $(date) — DONE"
            exit 0
            ;;
        *)
            echo "Usage: $0 [V3F|V5|V6]_[off|on]_[tp|fp] (single run) | no arg = full 12 runs"
            exit 1
            ;;
    esac
fi

for bam_tag in V3F V5 V6; do
    if [ ! -f "${BAMS[$bam_tag]}" ]; then
        echo "  [SKIP] $bam_tag BAM not available: ${BAMS[$bam_tag]}"
        continue
    fi
    echo
    echo "=== $bam_tag BAM ==="
    run_ism "$bam_tag" "off" "" "$VCF_TP" "tp"
    run_ism "$bam_tag" "off" "" "$VCF_FP" "fp"
    run_ism "$bam_tag" "on" "--germline-hp-only" "$VCF_TP" "tp"
    run_ism "$bam_tag" "on" "--germline-hp-only" "$VCF_FP" "fp"
done

echo
echo "[Phase C with significance] $(date) — DONE"
