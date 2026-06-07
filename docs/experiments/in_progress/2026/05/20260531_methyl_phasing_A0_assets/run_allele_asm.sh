#!/bin/bash
# 單 SNP allele-ASM AUC：tumor + normal 配對（同區域）。唯讀 BAM。
set -uo pipefail
cd "$(dirname "$0")"
CHRS="chr1 chr7 chr8 chr15 chr20 chr22"
echo "[allele-asm] start $(date '+%F %T') chrs=$CHRS"
# 組合 sample×chr 任務，-P 8 平行（tumor big7 / normal big8 不同碟）
for s in tumor normal; do for c in $CHRS; do echo "$s $c"; done; done | \
  xargs -P 8 -n 2 bash -c '
    s="$0"; c="$1"; t0=$(date +%s)
    timeout 2400 python3 allele_asm_auc.py --sample "$s" --chrom "$c" --max-sites 120 --shuffleK 100 \
      --out-dir . >/dev/null 2>>allele_asm_errors.log \
      && echo "[ok] $s $c $(( $(date +%s)-t0 ))s" || echo "[FAIL] $s $c $(( $(date +%s)-t0 ))s"
  '
echo "[allele-asm] ALL DONE $(date '+%F %T')"
ls -1 allele_asm_tumor_*.json allele_asm_normal_*.json 2>/dev/null | wc -l | xargs echo "JSON count:"
