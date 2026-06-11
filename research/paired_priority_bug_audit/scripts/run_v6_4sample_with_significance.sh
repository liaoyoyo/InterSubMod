#!/bin/bash
# Phase 2 Cycle 2 Step B1' — V6 4 samples ISM rerun WITH significance
#
# Purpose:
#   Existing phaseD V6 BAM ISM outputs have only header-row significance_summary.csv
#   (--no-distance-matrix was set in phaseD). For Cycle 2 cross-sample augmented
#   filter validation we need full significance (incl. distance / permanova /
#   methyl families). This script re-runs ISM for 4 samples × V6 BAM × off flag
#   × {tp,fp} = 8 runs, with significance enabled.
#
# Design choices (per Cycle 2 plan v2.0):
#   - V6 BAM only (Cycle 1 picked V6 over V3F/V5)
#   - off flag only (Cycle 1 Track A used V6_off_*)
#   - tp + fp labels for ΔF1 computation downstream
#   - 4 samples: H1437 / H2009 / HCC1954 / HCC1937
#   - Adds HCC1395 = 5 samples total in B3 cross-sample apply (HCC1395 reused
#     from existing phaseC_genome_three_way_with_significance V6_off_{tp,fp})
#
# Usage:
#   ./run_v6_4sample_with_significance.sh                      # full 8 runs
#   ./run_v6_4sample_with_significance.sh H1437 tp             # pilot single run
#   ./run_v6_4sample_with_significance.sh H1437 fp             # single FP
#
# Output:
#   $REPO/research/paired_priority_bug_audit/phaseC_v6_4sample_with_significance/
#       {sample}/V6_off_{tp,fp}/significance_summary.csv (+ distance matrices)
#
# Author: Phase 2 Cycle 2 Step B1' (2026-05-18)
# Plan: /bip7_disk/liaoyoyo2001/.claude/plans/v6-optimized-wadler.md

set -euo pipefail

# --- Disk safety (memory: feedback_tmp_disk_full_pipeline_pitfall.md) ---
export TMPDIR=/big7_disk/liaoyoyo2001/tmp
mkdir -p "$TMPDIR"

# --- Paths ---
REPO=/big7_disk/liaoyoyo2001/InterSubMod
ISM=$REPO/build/bin/inter_sub_mod
REF=/big7_disk/liaoyoyo2001/data/ref/GRCh38_no_alt_analysis_set.fasta
THREADS=24
WINDOW=5000

# Per-sample resources (Agent B1 deep search verified 2026-05-18)
declare -A BAMS
declare -A LOH_BEDS
BAMS[H1437]=/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/H1437/tumor_tagged.bam
BAMS[H2009]=/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/H2009/tumor_tagged.bam
BAMS[HCC1954]=/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/HCC1954/tumor_tagged.bam
BAMS[HCC1937]=/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/HCC1937/tumor_tagged.bam

LOH_BEDS[H1437]=/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/H1437/tumor_phased_LOH.bed
LOH_BEDS[H2009]=/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/H2009/tumor_phased_LOH.bed
LOH_BEDS[HCC1954]=/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/HCC1954/tumor_phased_LOH.bed
LOH_BEDS[HCC1937]=/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/HCC1937/tumor_phased_LOH.bed

# VCFs on big8 (read-only)
VCF_BASE=/big8_disk/liaoyoyo2001/Somatic_methylation_analysis_output

OUT_BASE=$REPO/research/paired_priority_bug_audit/phaseC_v6_4sample_with_significance
mkdir -p "$OUT_BASE"

# --- run_ism: single ISM invocation ---
run_ism() {
    local sample=$1
    local label=$2          # tp | fp
    local vcf_path=$VCF_BASE/$sample/VCF/filtered_snv_${label}.vcf.gz
    local out_dir=$OUT_BASE/$sample/V6_off_${label}
    mkdir -p "$out_dir"

    if [ ! -f "$vcf_path" ]; then
        echo "  [SKIP] VCF missing: $vcf_path"
        return 1
    fi
    if [ ! -f "${BAMS[$sample]}" ]; then
        echo "  [SKIP] BAM missing: ${BAMS[$sample]}"
        return 1
    fi
    if [ ! -f "${LOH_BEDS[$sample]}" ]; then
        echo "  [SKIP] LOH missing: ${LOH_BEDS[$sample]}"
        return 1
    fi

    local n_snvs=$(zcat "$vcf_path" | grep -v '^#' | wc -l)
    echo "  [run] sample=$sample label=$label n_snvs=$n_snvs → $out_dir"
    local t0=$(date +%s)

    $ISM \
        --tumor-bam "${BAMS[$sample]}" \
        --reference "$REF" \
        --vcf "$vcf_path" \
        --loh-bed "${LOH_BEDS[$sample]}" \
        --output-dir "$out_dir" \
        --window-size $WINDOW \
        --threads $THREADS \
        --log-level warn \
        > "$out_dir/run.log" 2>&1

    local rc=$?
    local t1=$(date +%s)
    local elapsed=$((t1-t0))
    if [ $rc -ne 0 ]; then
        echo "    [FAIL] rc=$rc elapsed=${elapsed}s — see $out_dir/run.log"
        return $rc
    fi

    # Quick sanity check
    local nrows=0
    if [ -f "$out_dir/significance_summary.csv" ]; then
        nrows=$(wc -l < "$out_dir/significance_summary.csv")
    fi
    echo "    [OK] elapsed=${elapsed}s significance_rows=$nrows"
}

# --- Argument parsing ---
SAMPLES=(H1437 H2009 HCC1954 HCC1937)
LABELS=(tp fp)

if [ $# -ge 1 ]; then
    if [ $# -ne 2 ]; then
        echo "Usage: $0 [sample] [tp|fp]"
        echo "       $0           # full 8 runs"
        exit 1
    fi
    SAMPLES=("$1")
    LABELS=("$2")
fi

echo "[Phase 2 Cycle 2 B1'] $(date)"
echo "  Samples: ${SAMPLES[*]}"
echo "  Labels:  ${LABELS[*]}"
echo "  OUT_BASE: $OUT_BASE"
echo "  TMPDIR:   $TMPDIR"
echo

# --- Main loop ---
fail_count=0
for sample in "${SAMPLES[@]}"; do
    echo "=== $sample ==="
    for label in "${LABELS[@]}"; do
        if ! run_ism "$sample" "$label"; then
            fail_count=$((fail_count+1))
        fi
    done
    echo
done

echo "[Phase 2 Cycle 2 B1'] $(date) — DONE  (failures: $fail_count)"
exit $fail_count
