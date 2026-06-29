#!/usr/bin/env bash
# per-sample clone/subclone 樹 pipeline driver(2026-06-29 多樣本任務)。
# 用 env-var 化的 sm_*.py 跑:sm_linkage(5-way)→merge→multilocus(5-way)→region_integration。
# 用法: run_sample_subclone.sh SAMPLE TBAM NBAM VD VCF_MODE CNBED WORKDIR
set -uo pipefail
SAMPLE="$1"; TBAM="$2"; NBAM="$3"; VD="$4"; VCF_MODE="$5"; CNBED="$6"; WORKDIR="$7"
PY=/bip7_disk/liaoyoyo2001/miniconda3/bin/python3
SCD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$WORKDIR"
export SM_TBAM="$TBAM" SM_NBAM="$NBAM" SM_VD="$VD" SM_VCF_MODE="$VCF_MODE" SM_CNBED="$CNBED" SM_WORKDIR="$WORKDIR" SM_DATA="$WORKDIR"
echo "[$SAMPLE] start $(date '+%H:%M:%S')  TBAM=$TBAM  VCF_MODE=$VCF_MODE  CNBED=${CNBED:-<none>}"

# Stage 1: sm_linkage 5-way parallel
declare -a SPLIT=("chr1,chr6,chr11,chr16,chr21:1" "chr2,chr7,chr12,chr17,chr22:2" "chr3,chr8,chr13,chr18:3" "chr4,chr9,chr14,chr19:4" "chr5,chr10,chr15,chr20:5")
for s in "${SPLIT[@]}"; do
  $PY "$SCD/sm_linkage_genomewide.py" "${s%:*}" "$WORKDIR/sm_linkage_gw_part_${s#*:}.json" > "$WORKDIR/linkage_part_${s#*:}.log" 2>&1 &
done
wait
echo "[$SAMPLE] linkage parts done $(date '+%H:%M:%S')"
ls "$WORKDIR"/sm_linkage_gw_part_*.json >/dev/null 2>&1 || { echo "[$SAMPLE] FAIL: no linkage parts"; exit 1; }

# Stage 2: merge
$PY "$SCD/sm_merge_parts.py" > "$WORKDIR/merge.log" 2>&1 || { echo "[$SAMPLE] FAIL merge"; cat "$WORKDIR/merge.log"; exit 1; }
echo "[$SAMPLE] merged $(date '+%H:%M:%S')"

# Stage 2b: pair_lists (產 lists/{config}_{hp}.tsv + per_sSNV_census.tsv;topology 需 HP/TP-FP 標註)
$PY "$SCD/sm_pair_lists.py" > "$WORKDIR/pair_lists.log" 2>&1 || { echo "[$SAMPLE] FAIL pair_lists"; cat "$WORKDIR/pair_lists.log"; exit 1; }
echo "[$SAMPLE] pair_lists done $(date '+%H:%M:%S')"

# Stage 3: multilocus 5-way parallel
for s in "${SPLIT[@]}"; do
  $PY "$SCD/sm_multilocus_combinations.py" "${s%:*}" "$WORKDIR/ml_part_${s#*:}.json" > "$WORKDIR/ml_part_${s#*:}.log" 2>&1 &
done
wait
echo "[$SAMPLE] multilocus done $(date '+%H:%M:%S')"
ls "$WORKDIR"/ml_part_*.json >/dev/null 2>&1 || { echo "[$SAMPLE] FAIL: no ml_part"; exit 1; }

# Stage 4: region integration
$PY "$SCD/sm_region_integration.py" > "$WORKDIR/region_integration.log" 2>&1 || { echo "[$SAMPLE] FAIL region_integration"; cat "$WORKDIR/region_integration.log"; exit 1; }
echo "[$SAMPLE] region_integration done $(date '+%H:%M:%S')"

# Stage 5: topology (SM_DATA=WORKDIR;讀 sm_region_integration.json + lists/;chr17 缺則 None)
$PY "$SCD/topology_analysis.py" > "$WORKDIR/topology.log" 2>&1 || { echo "[$SAMPLE] FAIL topology"; cat "$WORKDIR/topology.log"; exit 1; }
echo "[$SAMPLE] topology done $(date '+%H:%M:%S')"

# Stage 6: candidate_scoring
$PY "$SCD/candidate_scoring.py" > "$WORKDIR/candidate_scoring.log" 2>&1 || { echo "[$SAMPLE] FAIL candidate_scoring"; cat "$WORKDIR/candidate_scoring.log"; exit 1; }
echo "[$SAMPLE] DONE all $(date '+%H:%M:%S') -> $WORKDIR/{sm_region_integration,topology_per_region,candidate_scoring}.json"
ls -la "$WORKDIR"/topology_per_region.json "$WORKDIR"/candidate_scoring.json 2>/dev/null
