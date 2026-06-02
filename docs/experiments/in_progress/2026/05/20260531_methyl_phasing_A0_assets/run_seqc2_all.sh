#!/bin/bash
# 平行跑全 autosome seqc2_cn_methyl.py（產數字階段，確定性，唯讀 BAM）
set -uo pipefail
cd "$(dirname "$0")"
CHRS=$(seq 1 22 | sed 's/^/chr/')
echo "[seqc2-all] launching $(echo $CHRS | wc -w) chromosomes (parallel 6, max 30/status, K=50)..."
echo "$CHRS" | tr ' ' '\n' | xargs -P 6 -I {} python3 seqc2_cn_methyl.py --chrom {} --max-per-status 30 --K 50
echo "[seqc2-all] all chromosomes done."
ls -1 seqc2_cn_methyl_chr*.json 2>/dev/null | wc -l | xargs echo "JSON files produced:"
