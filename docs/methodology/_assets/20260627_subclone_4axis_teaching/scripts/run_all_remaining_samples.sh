#!/usr/bin/env bash
# 序列跑其餘 5 樣本(COLO829 已 pilot 完;HCC1395 凍結已有)→ clone/subclone 樹。
# 每樣本 ~45min,序列(非並行,避免大 BAM 撞 IO/mem,§8)。CN=unknown(SAVANA 後補)。
set -uo pipefail
C=/big7_disk/liaoyoyo2001/big7_disk_output/canonical
OUT=/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone
SCD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DRV="$SCD/run_sample_subclone.sh"
# sample : normal_bam
declare -A NB=(
 [H1437]=/big8_disk/data/H1437/ONT/H1437BL.bam
 [H2009]=/big8_disk/data/H2009/ONT/H2009BL.bam
 [HCC1395_DORADO]=/big8_disk/liaoyoyo2001/data/bam/HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam
 [HCC1937]=/big8_disk/data/HCC1937/ONT/HCC1937BL.bam
 [HCC1954]=/big8_disk/data/HCC1954/ONT/HCC1954BL.bam
)
for S in H1437 H2009 HCC1395_DORADO HCC1937 HCC1954; do
  LP=$(ls -d $C/$S/paired_full/*complete_matrix/longphase_s 2>/dev/null | head -1)
  TB="$LP/${S}_tagged.bam"
  if [ ! -f "$TB" ] || [ ! -f "${NB[$S]}" ]; then echo "[$S] SKIP: missing input (TB=$TB NB=${NB[$S]})"; continue; fi
  echo "===== [$S] $(date '+%F %H:%M:%S') ====="
  bash "$DRV" "$S" "$TB" "${NB[$S]}" "$LP" single "" "$OUT/$S" 2>&1
  echo "===== [$S] end $(date '+%H:%M:%S') ====="
done
echo "ALL REMAINING SAMPLES DONE $(date '+%F %H:%M:%S')"
