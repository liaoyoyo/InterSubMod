#!/bin/bash
set -uo pipefail
cd "$(dirname "$0")"
echo "[inv] start $(date '+%F %T')"
echo "chr15 chr20 chr22 chr19" | tr ' ' '\n' | xargs -P 4 -I {} bash -c '
  c="$1"; t0=$(date +%s)
  timeout 3000 python3 unphase_inventory.py --chrom "$c" --out-dir . 2>/dev/null \
    && echo "[ok] $c $(( $(date +%s)-t0 ))s" || echo "[FAIL] $c $(( $(date +%s)-t0 ))s"
' _ {}
echo "[inv] ALL DONE $(date '+%F %T')"
