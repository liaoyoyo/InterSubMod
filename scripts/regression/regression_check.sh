#!/bin/bash
# ISM regression check — 確認 C++ 改動沒改變結果（結果不變機制 / 主要研究基準守門）
# 固定 canonical 配置(chr1 ±5000 BERNOULLI paired) → 比 golden baseline → PASS/FAIL
# 用法:
#   regression_check.sh --update-golden   # 重建 golden(基準刻意改變/首次建立)
#   regression_check.sh                    # 比對(每次 C++ 改動後跑,確認結果不變)
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
BIN="$ROOT/build/bin/inter_sub_mod"
TUMOR=/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam
NORMAL=/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam
REF=/big8_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa
VCF=/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp_chr1.vcf.gz
GOLDEN="$HERE/golden_chr1_w5000_bernoulli.tsv"
OUT="$ROOT/output/_tmp_regression"
export TMPDIR=/big7_disk/liaoyoyo2001/tmp; mkdir -p "$TMPDIR"
rm -rf "$OUT"; mkdir -p "$OUT"

# 固定 canonical 配置（= 主要研究基準）
"$BIN" -t "$TUMOR" -n "$NORMAL" -r "$REF" -v "$VCF" -w 5000 -j 16 \
  --distance-metric BERNOULLI --no-output-distance-matrix -o "$OUT" > "$OUT/run.log" 2>&1

# 提取決定性結果欄（per region）
python3 -c "
import pandas as pd
d=pd.read_csv('$OUT/significance_summary.csv')
cols=['RegionID','NumReads','NumCpGs','GlobalP','CramersV','PassedGating','ClusterPermanovaValid','Significant','VerificationClass']
d[cols].sort_values('RegionID').to_csv('$OUT/current.tsv',sep='\t',index=False,float_format='%.6g')
"
if [ "${1:-}" = "--update-golden" ]; then
  cp "$OUT/current.tsv" "$GOLDEN"; echo "[regression] golden 已建/更新: $(($(wc -l < "$GOLDEN")-1)) 位點 @ commit $(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null)"
  rm -rf "$OUT"; exit 0
fi
[ ! -f "$GOLDEN" ] && { echo "[regression] 無 golden — 先跑 --update-golden"; exit 2; }
if diff -q "$GOLDEN" "$OUT/current.tsv" >/dev/null; then
  echo "[regression] ✅ PASS — 結果與 golden 完全一致 ($(($(wc -l < "$GOLDEN")-1)) 位點)"; rm -rf "$OUT"; exit 0
else
  echo "[regression] 🔴 FAIL — 結果改變！差異(golden<→current>):"; diff "$GOLDEN" "$OUT/current.tsv" | head -20
  echo "(完整輸出留 $OUT 供調查)"; exit 1
fi
