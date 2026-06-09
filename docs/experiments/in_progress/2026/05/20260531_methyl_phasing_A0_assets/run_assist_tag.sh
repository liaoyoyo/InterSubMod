#!/bin/bash
set -uo pipefail
cd "$(dirname "$0")"
CHRS="chr1 chr7 chr15 chr20 chr22"
echo "[assist] start $(date '+%F %T') chrs=$CHRS"
echo "$CHRS" | tr ' ' '\n' | xargs -P 5 -I {} bash -c '
  c="$1"; t0=$(date +%s)
  timeout 2400 python3 assist_tag_separable.py --chrom "$c" --max-blocks 50 --out-dir . 2>/dev/null \
    && echo "[ok] $c $(( $(date +%s)-t0 ))s" || echo "[FAIL] $c $(( $(date +%s)-t0 ))s"
' _ {}
echo "[assist] ALL DONE $(date '+%F %T')"
ls -1 assist_tag_chr*.json 2>/dev/null | wc -l | xargs echo "JSON:"
