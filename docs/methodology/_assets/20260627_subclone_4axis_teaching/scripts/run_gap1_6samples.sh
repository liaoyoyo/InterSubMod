#!/usr/bin/env bash
# gap#1 partial-read subcube 救回 — 擴到其他 6 樣本(2026-07-05)。
# 每樣本: 備份舊 ml_part → multilocus(gap#1 builder,含 subread_groups,讀 tagged BAM)
#   → region_integration(傳 subread_groups) → topology(E_subcube_recovered) → candidate_scoring。
# 6 樣本 cn=unknown → SM_CNBED="" → m-通道自動不可用(cn==unknown 短路,不誤讀 HCC1395 整數 CN)。
# multilocus 3-way 平行(IO-bound,不 thrash disk),downstream 逐樣本(快,無 BAM)。
set -uo pipefail
SCD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY=python3
MSROOT=/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone
C=/big7_disk/liaoyoyo2001/big7_disk_output/canonical
CHROMS=$($PY -c "print(','.join('chr'+str(i) for i in range(1,23)))")
SAMPLES=(COLO829 H1437 H2009 HCC1395_DORADO HCC1937 HCC1954)

echo "=== STAGE multilocus (gap#1, 3-way parallel) $(date '+%H:%M:%S') ==="
run_ml() {
  local s="$1"; local WD="$MSROOT/$s"
  local lp; lp=$(ls -d $C/$s/paired_full/*complete_matrix/longphase_s 2>/dev/null | head -1)
  local TBAM="$lp/${s}_tagged.bam"
  mkdir -p "$WD/ml_part_backup_pregap1"; mv "$WD"/ml_part_*.json "$WD/ml_part_backup_pregap1/" 2>/dev/null || true
  SM_WORKDIR="$WD" SM_TBAM="$TBAM" SM_VD="$lp" SM_VCF_MODE="single" SM_CNBED="" \
    $PY "$SCD/sm_multilocus_combinations.py" "$CHROMS" "$WD/ml_part_gap1.json" \
    > "$WD/gap1_multilocus.log" 2>&1 && echo "[$s] multilocus DONE $(date '+%H:%M:%S')" || echo "[$s] multilocus FAIL"
}
export -f run_ml; export MSROOT C CHROMS SCD PY
# 3-way 平行
printf '%s\n' "${SAMPLES[@]}" | xargs -P 3 -I {} bash -c 'run_ml "$@"' _ {}

echo "=== STAGE downstream (region_integration→topology→candidate_scoring) $(date '+%H:%M:%S') ==="
for s in "${SAMPLES[@]}"; do
  WD="$MSROOT/$s"
  [ -f "$WD/ml_part_gap1.json" ] || { echo "[$s] SKIP downstream(no ml_part_gap1)"; continue; }
  SM_WORKDIR="$WD" SM_DATA="$WD" SM_VCF_MODE="single" SM_CNBED="" $PY "$SCD/sm_region_integration.py" > "$WD/gap1_region_integration.log" 2>&1 \
    && SM_DATA="$WD" $PY "$SCD/topology_analysis.py" > "$WD/gap1_topology.log" 2>&1 \
    && SM_DATA="$WD" $PY "$SCD/candidate_scoring.py" > "$WD/gap1_candidate_scoring.log" 2>&1 \
    && SM_DATA="$WD" $PY "$SCD/enumerate_candidate_trees.py" > "$WD/gap1_enumerate.log" 2>&1 \
    && echo "[$s] downstream DONE $(date '+%H:%M:%S')" || echo "[$s] downstream FAIL"
done
echo "=== ALL DONE $(date '+%H:%M:%S') ==="
