#!/bin/bash
# v2 全基因組（22 autosome）平行，吃 ~80% 資源（48 核 → 38 平行）。
# 含 (a) CpG-SNP 排除 (b) 全基因組 (c) LOH 個案細節。depth 已 pysam in-process。
# 大染色體 het 多 → 每 chr 一個 process，外層 38 平行讓 22 條同時起跑。
set -uo pipefail
cd "$(dirname "$0")"
CHRS=$(seq 1 22 | sed 's/^/chr/')
echo "[v2-genome] launch 22 chr, parallel 22 (each 1 core, ~80% of 48), max 25/status K=30  $(date '+%H:%M:%S')"
echo "$CHRS" | tr ' ' '\n' | xargs -P 22 -I {} bash -c '
  c="$1"
  t0=$(date +%s)
  timeout 5400 python3 seqc2_cn_methyl_v2.py --chrom "$c" --max-per-status 25 --K 30 \
    && echo "[v2-genome] $c OK $(( $(date +%s)-t0 ))s" \
    || echo "[v2-genome] $c TIMEOUT/ERR $(( $(date +%s)-t0 ))s"
' _ {}
echo "[v2-genome] ALL DONE $(date '+%H:%M:%S')"
ls -1 seqc2_cn_methyl_v2_chr*.json 2>/dev/null | wc -l | xargs echo "v2 JSON files:"
