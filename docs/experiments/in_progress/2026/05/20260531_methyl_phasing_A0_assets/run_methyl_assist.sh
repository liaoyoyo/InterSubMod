#!/bin/bash
set -uo pipefail
cd "$(dirname "$0")"
CHRS="chr1 chr7 chr8 chr15 chr20 chr22"
echo "[assist-T2T3] start $(date '+%F %T') chrs=$CHRS"
echo "$CHRS" | tr ' ' '\n' | xargs -P 6 -I {} bash -c '
  c="$1"; t0=$(date +%s)
  timeout 3000 python3 methyl_assist_targets.py --chrom "$c" --max-sites 200 --out-dir . 2>/dev/null \
    && echo "[ok] $c $(( $(date +%s)-t0 ))s" || echo "[FAIL] $c $(( $(date +%s)-t0 ))s"
' _ {}
echo "[assist-T2T3] ALL DONE $(date '+%F %T')"
ls -1 methyl_assist_chr*.json 2>/dev/null | wc -l | xargs echo "JSON:"
