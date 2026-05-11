#!/bin/bash
# Phase C: HCC1395 全基因組 V3F vs V5 vs V6 three-way head-to-head
# 12 runs (3 BAM × 2 flag × 2 label), genome-wide ClairS-TO TP=30,490 + FP=4,842
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

OUT_BASE=$REPO/research/paired_priority_bug_audit/phaseC_genome_three_way
mkdir -p "$OUT_BASE"

echo "[Phase C 全基因組三向] $(date)"
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
    if [ -n "$flag_arg" ]; then
        $ISM \
            --tumor-bam "${BAMS[$bam_tag]}" \
            --reference "$REF" \
            --vcf "$vcf_path" \
            --loh-bed "$LOH_BED" \
            --output-dir "$out_dir" \
            --window-size 5000 \
            --threads $THREADS \
            --no-distance-matrix \
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
            --no-distance-matrix \
            --log-level warn \
            > "$out_dir/run.log" 2>&1
    fi
    local t1=$(date +%s)
    echo "    elapsed: $((t1-t0))s"
}

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
echo "[Phase C] $(date) — DONE"
