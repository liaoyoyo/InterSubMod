#!/bin/bash
set -uo pipefail
cd "$(dirname "$0")"
CHRS="chr7 chr15 chr20 chr22"
echo "[rigor] start $(date '+%F %T')"
( for c in $CHRS; do echo "t1 $c"; echo "t3 $c"; echo "t2 $c"; done ) | xargs -P 8 -n 2 bash -c '
  s="$0"; c="$1"; t0=$(date +%s)
  case "$s" in
    t1) timeout 3000 python3 rigor_t1.py --chrom "$c" --max-blocks 60 --out-dir . 2>/dev/null ;;
    t2) timeout 3000 python3 rigor_t2.py --chrom "$c" --max-sites 200 --out-dir . 2>/dev/null ;;
    t3) timeout 3000 python3 rigor_t3.py --chrom "$c" --max-sites 250 --out-dir . 2>/dev/null ;;
  esac && echo "[ok] $s $c $(( $(date +%s)-t0 ))s" || echo "[FAIL] $s $c $(( $(date +%s)-t0 ))s"
'
echo "[rigor] ALL DONE $(date '+%F %T')"
ls -1 rigor_t1_chr*.json rigor_t2_chr*.json rigor_t3_chr*.json 2>/dev/null | wc -l | xargs echo "JSON:"
