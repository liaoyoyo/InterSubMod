#!/bin/bash
# Baseline ISM add-on: HCC1395 full-genome baseline × 2 flag × 2 label = 4 runs
# Drops into existing phaseC_genome_three_way/ as baseline_{off,on}_{tp,fp}
# Allows step1 build script to consume baseline alongside V3F/V5/V6 (→ 4-way)
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

BASELINE_BAM=/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_tagged.bam

OUT_BASE=$REPO/research/paired_priority_bug_audit/phaseC_genome_three_way
mkdir -p "$OUT_BASE"

echo "[Baseline ISM add-on] $(date)"
echo "  baseline BAM: $BASELINE_BAM"
echo "  TP VCF: $VCF_TP ($(zcat $VCF_TP | grep -v '^#' | wc -l) SNVs)"
echo "  FP VCF: $VCF_FP ($(zcat $VCF_FP | grep -v '^#' | wc -l) SNVs)"
echo "  LOH BED: $LOH_BED (note: V5 LOH bed reused — non-circular since LOH.bed same across V3F/V5/V6)"

run_ism() {
    local flag_label=$1
    local flag_arg=$2
    local vcf_path=$3
    local label=$4
    local out_dir=$OUT_BASE/baseline_${flag_label}_${label}
    mkdir -p "$out_dir"
    echo "  [run] BAM=baseline flag=$flag_label label=$label → $out_dir"
    local t0=$(date +%s)
    if [ -n "$flag_arg" ]; then
        $ISM \
            --tumor-bam "$BASELINE_BAM" \
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
            --tumor-bam "$BASELINE_BAM" \
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

echo
echo "=== baseline BAM (4 runs) ==="
run_ism "off" "" "$VCF_TP" "tp"
run_ism "off" "" "$VCF_FP" "fp"
run_ism "on" "--germline-hp-only" "$VCF_TP" "tp"
run_ism "on" "--germline-hp-only" "$VCF_FP" "fp"

echo
echo "[Baseline ISM add-on] $(date) — DONE"
