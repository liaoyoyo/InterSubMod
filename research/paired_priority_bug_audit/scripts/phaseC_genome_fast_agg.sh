#!/bin/bash
# Fast genome-wide aggregation: sum hp tag counts across all regions per BAM.
# Bypasses Python file-by-file overhead; parallel across runs.
set -euo pipefail

BASE=/big7_disk/liaoyoyo2001/InterSubMod/research/paired_priority_bug_audit/phaseC_genome_three_way

agg_run() {
    local run_dir=$1
    local label=$(basename "$run_dir")
    # Sum hp counts via cat | awk (fast: 1 fork per file but all reads concatenated)
    find "$run_dir" -name reads.tsv -print0 \
        | xargs -0 -P 4 -I{} awk -F'\t' 'NR>1{print $7}' {} \
        | sort | uniq -c | sort -rn > "${run_dir}/hp_summary.txt"
    echo "  ${label} done"
}

export -f agg_run

echo "[fast-agg] $(date) — aggregating hp counts per run"
for d in V3F V5 V6; do
    for f in off_tp off_fp on_tp on_fp; do
        agg_run "$BASE/${d}_${f}" &
    done
done
wait

echo "[fast-agg] $(date) — DONE"
echo
echo "=== Aggregated hp distribution (all 12 runs) ==="
for d in V3F V5 V6; do
    echo
    echo "--- ${d} ---"
    for f in off_tp off_fp on_tp on_fp; do
        echo "  ${f}:"
        cat "$BASE/${d}_${f}/hp_summary.txt" | sed 's/^/    /'
    done
done
