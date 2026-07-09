#!/usr/bin/env bash
# 全 7 樣本新骨幹重跑(2026-07-09 使用者定案:骨幹=ClairS PASS/LongPhase-S _sc.vcf,移除 is_somatic 粗重檢)。
# 每樣本:mlhp multilocus(5-way chrom 平行,SM_SOMATIC_VCF=somatic_pass)→ layered → region-view。
# 資源:每樣本 5-way + 2 樣本併發(~10 併發 BAM reader,平衡 disk I/O)。
# 評估:每階段 /usr/bin/time -v 記 wall time + peak RSS → profile.log。
set -uo pipefail
SCD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY=/bip7_disk/liaoyoyo2001/miniconda3/bin/python3
C=/big7_disk/liaoyoyo2001/big7_disk_output/canonical
MSROOT=/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone
PILOT=/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot
INTG=/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/sm_region_integration.json
PROF="${SCD}/../../20260618_subcluster_pilot/newbb_profile.log"
: > "$PROF"
declare -a SPLIT=("chr1,chr6,chr11,chr16,chr21:1" "chr2,chr7,chr12,chr17,chr22:2" "chr3,chr8,chr13,chr18:3" "chr4,chr9,chr14,chr19:4" "chr5,chr10,chr15,chr20:5")
SEQC2_GAIN=/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnv_gain_cn.bed
SEQC2_LOSS=/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnv_loss_cn.bed
SEQC2_CNBED=/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed
TIME=/usr/bin/time

prof(){ echo "[$(date '+%H:%M:%S')] $*" >> "$PROF"; }

run_sample(){
  local s="$1"; local WD CNBED CNG CNL CENSUS
  local -a SPLIT=("chr1,chr6,chr11,chr16,chr21:1" "chr2,chr7,chr12,chr17,chr22:2" "chr3,chr8,chr13,chr18:3" "chr4,chr9,chr14,chr19:4" "chr5,chr10,chr15,chr20:5")
  if [ "$s" = "HCC1395" ]; then WD="$PILOT"; CNBED="$SEQC2_CNBED"; CNG="$SEQC2_GAIN"; CNL="$SEQC2_LOSS";
  else WD="$MSROOT/$s"; CNBED=""; CNG=""; CNL=""; fi   # 6 樣本 cn=unknown → L2 CN 不可用
  local lp; lp=$(ls -d $C/$s/paired_full/*complete_matrix/longphase_s 2>/dev/null | head -1)
  local TBAM="$lp/${s}_tagged.bam"; local SPASS="$lp/somatic_pass.vcf.gz"
  [ -f "$TBAM" ] && [ -f "$SPASS" ] || { echo "[$s] FAIL: 缺 TBAM/somatic_pass"; prof "[$s] FAIL no input"; return; }
  local t0=$(date +%s)
  # backup 舊 mlhp
  mkdir -p "$WD/mlhp_oldbb_backup"; mv "$WD"/mlhp_part*.json "$WD/mlhp_oldbb_backup/" 2>/dev/null || true
  # Stage1 mlhp 5-way(新骨幹)
  for sp in "${SPLIT[@]}"; do
    SM_WORKDIR="$WD" SM_TBAM="$TBAM" SM_VD="$lp" SM_SOMATIC_VCF="$SPASS" SM_CNBED="$CNBED" \
      $TIME -v $PY "$SCD/sm_multilocus_combinations.py" "${sp%:*}" "$WD/mlhp_part_${sp#*:}.json" \
      > "$WD/mlhp_nb_${sp#*:}.log" 2>&1 &
  done
  wait
  ls "$WD"/mlhp_part_*.json >/dev/null 2>&1 || { echo "[$s] FAIL mlhp"; prof "[$s] FAIL mlhp"; return; }
  # peak RSS(取 5 片最大) + wall(整 Stage1)
  local rss=$(grep -h "Maximum resident" "$WD"/mlhp_nb_*.log 2>/dev/null | grep -oE '[0-9]+' | sort -rn | head -1)
  local t1=$(date +%s)
  # Stage2 layered + region-view
  SM_ML_GLOB="$WD/mlhp_part_*.json" SM_OUT="$WD/layered_reconstruction_${s}.json" SM_VERIFY_EVERY=20 \
    SM_CN_INT_GAIN="$CNG" SM_CN_INT_LOSS="$CNL" \
    $TIME -v $PY "$SCD/layered_tree_reconstruction.py" > "$WD/layered_nb.log" 2>&1
  local lrss=$(grep "Maximum resident" "$WD/layered_nb.log" | grep -oE '[0-9]+' | head -1)
  local integ_arg=""; [ "$s" = "HCC1395" ] && integ_arg="$INTG"
  SM_LAYERED="$WD/layered_reconstruction_${s}.json" SM_OUT="$WD/layered_region_view_${s}.json" SM_SAMPLE="$s" \
    SM_SOMATIC_VCF="$SPASS" SM_INTEGRATION="$integ_arg" \
    $PY "$SCD/build_region_view.py" > "$WD/region_view_nb.log" 2>&1
  local t2=$(date +%s)
  prof "[$s] mlhp: $((t1-t0))s (peak RSS $((rss/1024))MB) | layered+rv: $((t2-t1))s (layered peak $((lrss/1024))MB) | total $((t2-t0))s"
  echo "[$s] DONE total $((t2-t0))s (mlhp $((t1-t0))s / down $((t2-t1))s) $(date '+%H:%M:%S')"
}
export -f run_sample prof; export SCD PY C MSROOT PILOT INTG PROF TIME SEQC2_GAIN SEQC2_LOSS SEQC2_CNBED
export SPLIT  # 注:array 不能直接 export;下方改用列印傳入

echo "=== 全 7 樣本新骨幹重跑 START $(date '+%H:%M:%S') ==="
prof "=== START $(date '+%F %H:%M:%S') 機器 $(nproc) cores ==="
# 2 樣本併發(每樣本內 5-way → ~10 併發 reader)
printf '%s\n' HCC1395 COLO829 H1437 H2009 HCC1395_DORADO HCC1937 HCC1954 | xargs -P 2 -I {} bash -c 'run_sample "$@"' _ {}
echo "=== ALL_7_NEWBB_DONE $(date '+%H:%M:%S') ==="
prof "=== ALL_7_NEWBB_DONE $(date '+%F %H:%M:%S') ==="
cat "$PROF"
