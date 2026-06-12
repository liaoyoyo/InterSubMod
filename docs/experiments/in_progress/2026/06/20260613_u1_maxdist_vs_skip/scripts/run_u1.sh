#!/bin/bash
# U1: MAX_DIST vs SKIP 全量對比（BERNOULLI ±5000 對齊 canonical，純彙總）
set -uo pipefail
export TMPDIR=/big7_disk/liaoyoyo2001/tmp; mkdir -p "$TMPDIR"
BIN=/big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod
TUMOR=/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam
NORMAL=/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam
REF=/big8_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa
VCFDIR=/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup
OUT=/big7_disk/liaoyoyo2001/InterSubMod/output/_tmp_u1_maxdist_vs_skip
G0=$(date +%s)
for strat in MAX_DIST SKIP; do
  for cls in tp fp; do
    echo "=== $strat $cls start $(date +%T) ==="
    t0=$(date +%s)
    "$BIN" -t "$TUMOR" -n "$NORMAL" -r "$REF" -v "$VCFDIR/filtered_snv_${cls}.vcf.gz" \
      -o "$OUT/${strat}_${cls}" -w 5000 -j 40 --distance-metric BERNOULLI \
      --nan-distance-strategy "$strat" --no-output-distance-matrix > "$OUT/${strat}_${cls}.log" 2>&1
    echo "=== $strat $cls done exit=$? $(($(date +%s)-t0))s ==="
  done
done
echo "U1 ALL DONE total $(($(date +%s)-G0))s"
