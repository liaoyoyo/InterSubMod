#!/bin/bash
set -uo pipefail
cd "$(dirname "$0")"
echo "[t3local] start $(date '+%F %T')"
echo "chr7 chr15 chr20 chr22" | tr ' ' '\n' | xargs -P 4 -I {} bash -c '
  c="$1"; t0=$(date +%s)
  timeout 3000 python3 t3_local_allele.py --chrom "$c" --max-sites 250 --out-dir . 2>/dev/null \
    && echo "[ok] $c $(( $(date +%s)-t0 ))s" || echo "[FAIL] $c $(( $(date +%s)-t0 ))s"
' _ {}
echo "[t3local] ALL DONE $(date '+%F %T')"
