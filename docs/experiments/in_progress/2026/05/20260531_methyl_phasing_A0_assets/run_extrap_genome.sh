#!/bin/bash
# 外推驗證全 autosome（PS-block held-out 救援模擬）。8 平行（每條快 ~30-60s，不像 v2 慢）。
set -uo pipefail
cd "$(dirname "$0")"
CHRS=$(seq 1 22 | sed 's/^/chr/')
echo "[extrap-genome] 22 chr, parallel 12, 8 blocks × 3 窗 K=30  $(date '+%H:%M:%S')"
echo "$CHRS" | tr ' ' '\n' | xargs -P 12 -I {} bash -c '
  c="$1"; t0=$(date +%s)
  timeout 1200 python3 extrapolation_validation.py --chrom "$c" --n-blocks 8 --windows-per-block 3 --K 30 2>/dev/null \
    && echo "[extrap] $c OK $(( $(date +%s)-t0 ))s" || echo "[extrap] $c TIMEOUT $(( $(date +%s)-t0 ))s"
' _ {}
echo "[extrap-genome] ALL DONE $(date '+%H:%M:%S')"
ls -1 extrapolation_chr*.json 2>/dev/null | wc -l | xargs echo "extrap JSON:"
