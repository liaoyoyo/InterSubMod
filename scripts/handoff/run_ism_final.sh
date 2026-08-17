#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-only
#
# 對最終交付的 BAM 逐染色體跑 ISM 甲基分析。
#
#   setsid nohup bash run_ism_final.sh > <LOG> 2>&1 &
#
# 參數為什麼是這幾個（皆為刻意選擇，不是預設值照抄）：
#   -w 5000               同-BAM 配對 A/B 實測：每對共同 CpG 20.35→59.80（2.94×）、
#                         MultiGroupNoLabel −57%、Noise −58%。窄窗（±1000）雖然
#                         valid-pair ratio 較高，但那是分母效應（納入的 read 較少）。
#   --distance-metric BERNOULLI
#                         Config.hpp:40 的程式預設，也是 regression golden
#                         (golden_chr1_w5000_bernoulli_skip.tsv) 的基準。
#                         先前某輪用 NHD 是未記錄的偏離。
#   --nan-distance-strategy SKIP
#                         MAX_DIST 會用 1.0 假距離污染矩陣，系統性低估甲基分群。
#   --group-by-tag HP,ALT,lc,lu,lv
#                         五個軸各自獨立檢定，讓每個位點能回報「是哪個軸解釋了
#                         甲基差異」，而不是只給一個混合的 p 值。
#   --require-tag-status U
#                         只採唯一解的 read。放寬到含 '+' 的部分覆蓋 read 會引入
#                         長度 confound：實測帶 '+' 的 read 中位 18,397 bp vs
#                         不帶 '+' 的 28,421 bp（1.54×），而 read 長度與甲基有關聯。

set -uo pipefail

BIN=${BIN:-/big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod}
FINAL=${FINAL:-/bip7_disk/liaoyoyo2001/HCC1395_final}
SAMPLE=${SAMPLE:-HCC1395}
NORMAL=${NORMAL:-/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam}
REF=${REF:-/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta}

BAM="$FINAL/bam/${SAMPLE}.lineage_tagged.bam"
VCFD="$FINAL/inputs/vcf"
OUT="$FINAL/ism"
mkdir -p "$OUT/logs"

[[ -f "$BAM" ]] || { echo "缺 BAM：$BAM"; exit 2; }
[[ -d "$VCFD" ]] || { echo "缺 VCF 目錄：$VCFD"; exit 2; }

echo "=== ISM $(date -Is) ==="
echo "  tumor  $BAM"
echo "  normal $NORMAL"
echo "  vcf    $VCFD"
fail=0
for c in chr{1..22}; do
    [[ -f "$VCFD/${c}.vcf.gz" ]] || { echo "  $c SKIP（無 VCF）"; continue; }
    # 已完成的跳過 —— 讓這支也可續跑（同 retag 的教訓）
    if [[ -f "$OUT/$c/significance_summary.csv" && -f "$OUT/$c/run_summary.json" ]]; then
        echo "  $c SKIP（已完成）"; continue
    fi
    t0=$SECONDS
    if "$BIN" -t "$BAM" -n "$NORMAL" -r "$REF" -v "$VCFD/${c}.vcf.gz" -o "$OUT/$c" \
        -w 5000 -j 24 \
        --distance-metric BERNOULLI \
        --min-common-coverage 3 --nan-distance-strategy SKIP \
        --compute-distance-matrix --output-distance-matrix \
        --group-by-tag HP,ALT,lc,lu,lv --require-tag-status U \
        > "$OUT/logs/${c}.log" 2>&1; then
        n=$(tail -n +2 "$OUT/$c/significance_summary.csv" 2>/dev/null | wc -l)
        echo "  $c OK  region=$n  ($((SECONDS-t0))s)"
    else
        echo "  $c FAIL  見 $OUT/logs/${c}.log"; fail=$((fail+1))
    fi
done
tot=$(cat "$OUT"/chr*/significance_summary.csv 2>/dev/null | grep -vc '^RegionID' || true)
echo "=== 完成 $(date -Is)　失敗 $fail／22　region 合計 $tot ==="
