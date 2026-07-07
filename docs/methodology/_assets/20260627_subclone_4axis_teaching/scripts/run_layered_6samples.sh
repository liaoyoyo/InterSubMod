#!/usr/bin/env bash
# 5b:6 樣本新分層 pipeline(2026-07-07)。mlhp multilocus(per-HP-家族,3-way平行讀 tagged BAM)
#   → layered driver → region-view。config 同 run_gap1_6samples(SM_VCF_MODE=single/SM_CNBED=""/cn=unknown)。
# 6 樣本 cn=unknown → L2 CN 通道不可用(誠實);L0/L1/L3 可跑。
set -uo pipefail
SCD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY=python3
MSROOT=/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone
C=/big7_disk/liaoyoyo2001/big7_disk_output/canonical
CHROMS=$($PY -c "print(','.join('chr'+str(i) for i in range(1,23)))")
SAMPLES=(COLO829 H1437 H2009 HCC1395_DORADO HCC1937 HCC1954)

echo "=== STAGE mlhp multilocus (per-HP-家族, 3-way 平行) $(date '+%H:%M:%S') ==="
run_ml() {
  local s="$1"; local WD="$MSROOT/$s"
  local lp; lp=$(ls -d $C/$s/paired_full/*complete_matrix/longphase_s 2>/dev/null | head -1)
  local TBAM="$lp/${s}_tagged.bam"
  [ -f "$TBAM" ] || { echo "[$s] mlhp FAIL: no tagged BAM ($TBAM)"; return; }
  [ -f "$WD/sm_linkage_genomewide.json" ] || { echo "[$s] mlhp FAIL: no census"; return; }
  SM_WORKDIR="$WD" SM_TBAM="$TBAM" SM_VD="$lp" SM_VCF_MODE="single" SM_CNBED="" \
    $PY "$SCD/sm_multilocus_combinations.py" "$CHROMS" "$WD/mlhp_part.json" \
    > "$WD/mlhp_multilocus.log" 2>&1 && echo "[$s] mlhp DONE $(date '+%H:%M:%S')" || echo "[$s] mlhp FAIL"
}
export -f run_ml; export MSROOT C CHROMS SCD PY
printf '%s\n' "${SAMPLES[@]}" | xargs -P 3 -I {} bash -c 'run_ml "$@"' _ {}

echo "=== STAGE layered + region-view (逐樣本) $(date '+%H:%M:%S') ==="
for s in "${SAMPLES[@]}"; do
  WD="$MSROOT/$s"
  [ -f "$WD/mlhp_part.json" ] || { echo "[$s] SKIP downstream(no mlhp_part)"; continue; }
  SM_ML="$WD/mlhp_part.json" SM_OUT="$WD/layered_reconstruction_${s}.json" SM_VERIFY_EVERY=20 \
    $PY "$SCD/layered_tree_reconstruction.py" > "$WD/layered.log" 2>&1 \
  && SM_LAYERED="$WD/layered_reconstruction_${s}.json" SM_OUT="$WD/layered_region_view_${s}.json" SM_SAMPLE="$s" \
     SM_CENSUS="$WD/sm_linkage_genomewide.json" \
     $PY "$SCD/build_region_view.py" > "$WD/region_view.log" 2>&1 \
  && echo "[$s] layered+regionview DONE $(date '+%H:%M:%S')" || echo "[$s] downstream FAIL(see $WD/layered.log)"
done
echo "=== ALL_6_DONE $(date '+%H:%M:%S') ==="
