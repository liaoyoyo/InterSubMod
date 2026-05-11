#!/bin/bash
# V6-C Phase B (V6 BAM): chr19 子集 ISM × 2 flag (germline-hp-only on/off)
# Three-way head-to-head with V3F + V5 BAMs (same protocol)
set -euo pipefail

export TMPDIR=/big7_disk/liaoyoyo2001/tmp
mkdir -p $TMPDIR

REPO=/big7_disk/liaoyoyo2001/InterSubMod
ISM=$REPO/build/bin/inter_sub_mod
TUMOR_BAM=/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_germline_absent_revert/tumor_tagged.bam
REF=/big7_disk/liaoyoyo2001/data/ref/GRCh38_no_alt_analysis_set.fasta
LOH_BED=/big7_disk/liaoyoyo2001/longphase-to-mod/output/threshold_compare/v5_flag/tumor_phased_LOH.bed
VCF_TP=/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp_chr19.vcf.gz
VCF_FP=/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_fp_chr19.vcf.gz
OUT_BASE=$REPO/research/paired_priority_bug_audit/v6c_phaseB_v6bam_runs
THREADS=24

mkdir -p "$OUT_BASE"
echo "[V6-C Phase B V6 BAM] $(date) — chr19 ISM × 2 flag (V6 BAM, three-way head-to-head)"
echo "  TUMOR_BAM (V6):  $TUMOR_BAM"
echo "  THREADS:         $THREADS"

run_ism() {
    local flag_label=$1
    local flag_arg=$2
    local vcf_path=$3
    local label=$4
    local out_dir=$OUT_BASE/${flag_label}_${label}
    mkdir -p "$out_dir"
    echo "  [run] flag=$flag_label label=$label → $out_dir"
    local t0=$(date +%s)
    if [ -n "$flag_arg" ]; then
        $ISM \
            --tumor-bam "$TUMOR_BAM" \
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
            --tumor-bam "$TUMOR_BAM" \
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
echo "=== Run 1/4: V6 BAM flag=off, TP ==="
run_ism "off" "" "$VCF_TP" "tp"

echo
echo "=== Run 2/4: V6 BAM flag=off, FP ==="
run_ism "off" "" "$VCF_FP" "fp"

echo
echo "=== Run 3/4: V6 BAM flag=on, TP ==="
run_ism "on" "--germline-hp-only" "$VCF_TP" "tp"

echo
echo "=== Run 4/4: V6 BAM flag=on, FP ==="
run_ism "on" "--germline-hp-only" "$VCF_FP" "fp"

echo
echo "[V6-C Phase B V6 BAM] $(date) — DONE"
