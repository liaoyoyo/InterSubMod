#!/bin/bash
set -uo pipefail
cd "$(dirname "$0")"
echo "chr7 chr15 chr20 chr22" | tr ' ' '\n' | xargs -P 4 -I {} bash -c '
  c="$1"; timeout 2400 python3 t3_reconcile.py --chrom "$c" --max-sites 200 --out-dir . 2>/dev/null && echo "[ok] $c" || echo "[FAIL] $c"' _ {}
echo DONE
