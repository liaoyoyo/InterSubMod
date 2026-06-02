#!/bin/bash
# 補充 loss-rich 染色體（循序，避免大染色體平行卡死）。已跑 chr5/6/7/21 跳過。
# depth 已改 pysam in-process；max-per-status 12 控時間。
set -uo pipefail
cd "$(dirname "$0")"
# loss 多的染色體（排除已完成）：chr2 chr18 chr11 chr13 + chr8 chr1（補 gain/loh/neutral 跨染色體一致性）
for c in chr2 chr18 chr11 chr13 chr8; do
  echo "[supp] start $c $(date '+%H:%M:%S')"
  timeout 1800 python3 seqc2_cn_methyl.py --chrom "$c" --max-per-status 12 --K 30
  echo "[supp] done $c $(date '+%H:%M:%S')"
done
echo "[supp] ALL DONE"
ls -1 seqc2_cn_methyl_chr*.json | wc -l | xargs echo "total JSON:"
