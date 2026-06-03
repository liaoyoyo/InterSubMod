#!/bin/bash
# 全 LOH 區 tumor VAF 系統分布（多染色體平行）。唯讀 BAM。
set -uo pipefail
cd "$(dirname "$0")"
# 有 LOH 的染色體（含已分析的 chr8/15/22 + 大染色體）
CHRS="chr1 chr2 chr5 chr8 chr11 chr13 chr15 chr16 chr18 chr20 chr22"
echo "[lohvaf-genome] $(echo $CHRS|wc -w) chr, parallel 11, max 300 sites/chr  $(date '+%H:%M:%S')"
echo "$CHRS" | tr ' ' '\n' | xargs -P 11 -I {} bash -c '
  c="$1"; t0=$(date +%s)
  timeout 1800 python3 loh_tumor_vaf_systematic.py --chrom "$c" --max-sites 300 2>/dev/null \
    && echo "[lohvaf] $c OK $(( $(date +%s)-t0 ))s" || echo "[lohvaf] $c TIMEOUT $(( $(date +%s)-t0 ))s"
' _ {}
echo "[lohvaf-genome] ALL DONE $(date '+%H:%M:%S')"
ls -1 loh_tumor_vaf_chr*.json 2>/dev/null | wc -l | xargs echo "VAF JSON:"
