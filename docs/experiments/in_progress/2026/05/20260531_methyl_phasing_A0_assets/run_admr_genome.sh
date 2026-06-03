#!/bin/bash
set -uo pipefail; cd "$(dirname "$0")"
CHRS="chr1 chr2 chr5 chr7 chr8 chr11 chr13 chr15 chr16 chr19 chr20 chr22"
echo "[admr-genome] $(echo $CHRS|wc -w) chr parallel 12, 100窗/chr $(date '+%H:%M:%S')"
echo "$CHRS"|tr ' ' '\n'|xargs -P 12 -I {} bash -c '
  c="$1";t0=$(date +%s)
  timeout 1500 python3 admr_loh_enrichment.py --chrom "$c" --n-windows 100 2>/dev/null \
    && echo "[admr] $c OK $(( $(date +%s)-t0 ))s" || echo "[admr] $c TIMEOUT"
' _ {}
echo "[admr-genome] DONE $(date '+%H:%M:%S')"
