#!/bin/bash
# 全量取樣門檻觀察 — 重跑 HCC1395 TP+FP+FN 開矩陣輸出（框架實例 #1 的 compute 步驟）
# 輸出 per-region 矩陣到臨時區（分析完刪）。用法: bash run_sampling_obs.sh
set -uo pipefail

export TMPDIR=/big7_disk/liaoyoyo2001/tmp   # 避免寫爆 root /tmp（feedback 教訓）
mkdir -p "$TMPDIR"

BIN=/big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod
TUMOR=/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam
NORMAL=/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam
REF=/big8_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa
VCFDIR=/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup
OUT=/big7_disk/liaoyoyo2001/InterSubMod/output/_tmp_sampling_obs_20260612_2201
THREADS=40

mkdir -p "$OUT"
GSTART=$(date +%s)
for cls in tp fp fn; do
  echo "=== [$cls] start $(date +%T) ==="
  cstart=$(date +%s)
  "$BIN" -t "$TUMOR" -n "$NORMAL" -r "$REF" \
    -v "$VCFDIR/filtered_snv_${cls}.vcf.gz" \
    -o "$OUT/$cls" -w 1000 -j "$THREADS" \
    --compute-distance-matrix --output-distance-matrix \
    --min-common-coverage 3 \
    > "$OUT/${cls}.run.log" 2>&1
  rc=$?
  cend=$(date +%s)
  ndir=$(ls -d "$OUT/$cls"/*/ 2>/dev/null | wc -l)
  echo "=== [$cls] done $(date +%T) | exit=$rc | $((cend-cstart))s | region dirs: $ndir ==="
done
GEND=$(date +%s)
echo "ALL DONE | total $((GEND-GSTART))s"
du -sh "$OUT" 2>/dev/null
