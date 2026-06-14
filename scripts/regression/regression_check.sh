#!/bin/bash
# ISM regression check — 結果不變機制（雙守護，2026-06-14 預設改 SKIP 後）
# 固定 canonical 配置(chr1 ±5000 BERNOULLI paired) → 比 golden → PASS/FAIL
#   SKIP(預設)     比 golden_chr1_w5000_bernoulli_skip.tsv   ← 主研究基準
#   MAX_DIST(對照) 比 golden_chr1_w5000_bernoulli.tsv        ← 確認改預設沒動 MAX_DIST 路徑(歷史名)
# 用法:
#   regression_check.sh --update-golden   # 重建兩個 golden
#   regression_check.sh                    # 雙路徑比對(每次 C++ 改動後跑)
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
BIN="$ROOT/build/bin/inter_sub_mod"
TUMOR=/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam
NORMAL=/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam
REF=/big8_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa
VCF=/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp_chr1.vcf.gz
GOLDEN_SKIP="$HERE/golden_chr1_w5000_bernoulli_skip.tsv"
GOLDEN_MAXDIST="$HERE/golden_chr1_w5000_bernoulli.tsv"
export TMPDIR=/big7_disk/liaoyoyo2001/tmp; mkdir -p "$TMPDIR"

run_extract() {  # $1=nan_strategy  $2=outdir
  rm -rf "$2"; mkdir -p "$2"
  "$BIN" -t "$TUMOR" -n "$NORMAL" -r "$REF" -v "$VCF" -w 5000 -j 16 \
    --distance-metric BERNOULLI --nan-distance-strategy "$1" --no-output-distance-matrix -o "$2" > "$2/run.log" 2>&1
  python3 -c "
import pandas as pd
d=pd.read_csv('$2/significance_summary.csv')
cols=['RegionID','NumReads','NumCpGs','GlobalP','CramersV','PassedGating','ClusterPermanovaValid','Significant','VerificationClass']
d[cols].sort_values('RegionID').to_csv('$2/current.tsv',sep='\t',index=False,float_format='%.6g')
"
}

OUT_SKIP="$ROOT/output/_tmp_regression_skip"
OUT_MAXD="$ROOT/output/_tmp_regression_maxdist"

if [ "${1:-}" = "--update-golden" ]; then
  run_extract SKIP "$OUT_SKIP";       cp "$OUT_SKIP/current.tsv" "$GOLDEN_SKIP"
  run_extract MAX_DIST "$OUT_MAXD";   cp "$OUT_MAXD/current.tsv" "$GOLDEN_MAXDIST"
  echo "[regression] golden 已建/更新: SKIP=$(($(wc -l < "$GOLDEN_SKIP")-1)) + MAX_DIST=$(($(wc -l < "$GOLDEN_MAXDIST")-1)) 位點 @ commit $(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null)"
  rm -rf "$OUT_SKIP" "$OUT_MAXD"; exit 0
fi

FAIL=0
check_one() {  # $1=label $2=strategy $3=golden $4=outdir
  [ ! -f "$3" ] && { echo "[regression] $1 無 golden ($3) — 先跑 --update-golden"; exit 2; }
  run_extract "$2" "$4"
  if diff -q "$3" "$4/current.tsv" >/dev/null; then
    echo "[regression] ✅ $1 PASS — 與 golden 一致 ($(($(wc -l < "$3")-1)) 位點)"; rm -rf "$4"
  else
    echo "[regression] 🔴 $1 FAIL — 結果改變！差異(golden<→current>):"; diff "$3" "$4/current.tsv" | head -15
    echo "(完整輸出留 $4)"; FAIL=1
  fi
}

check_one "SKIP(預設/主基準)" SKIP "$GOLDEN_SKIP" "$OUT_SKIP"
check_one "MAX_DIST(對照)" MAX_DIST "$GOLDEN_MAXDIST" "$OUT_MAXD"

if [ "$FAIL" -eq 0 ]; then
  echo "[regression] ✅ 雙守護 PASS (SKIP 預設 + MAX_DIST 對照)"; exit 0
else
  echo "[regression] 🔴 雙守護 FAIL"; exit 1
fi
