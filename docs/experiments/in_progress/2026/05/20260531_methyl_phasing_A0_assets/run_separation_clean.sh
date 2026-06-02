#!/bin/bash
# 乾淨序列跑：extract + separation，多個 anchor-充足 region。結果全寫檔。
set -uo pipefail
cd "$(dirname "$0")"
BAM=/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam
OUT=separation_results.txt
: > "$OUT"

run_region() {
  local name=$1 region=$2 mincpg=$3 mincol=$4
  echo "======== $name ($region) ========" >> "$OUT"
  python3 extract_per_read_methyl.py --bam "$BAM" --region "$region" --out-prefix "$name" --min-cpg "$mincpg" >> "$OUT" 2>&1
  python3 heatmap_and_separation.py --prefix "$name" --title "$name $region" --min-col-cov "$mincol" 2>> "$OUT" \
    | python3 -c "
import sys,json
d=json.load(sys.stdin); s=d.get('separation',{})
print('  RESULT reads=%s cpg=%s AUC=%s sil=%s null_p95=%s exceeds=%s verdict=%s' % (
  d.get('n_reads'), d.get('n_cpg_cols'),
  round(s['anchor_AUC_HP1_vs_HP2'],4) if s.get('anchor_AUC_HP1_vs_HP2') is not None else None,
  round(s['silhouette_HP1_vs_HP2'],4) if s.get('silhouette_HP1_vs_HP2') is not None else None,
  round(s['shuffle_null_silhouette_p95'],5) if s.get('shuffle_null_silhouette_p95') is not None else None,
  s.get('silhouette_exceeds_null_p95'), d.get('verdict')))
print('  hp_counts=%s' % d.get('hp_counts'))
" >> "$OUT" 2>&1
}

# anchor-充足 region：核心小窗（CpG 不過度分散），min-col-cov 放寬到 0.15
run_region brca2_core  "chr13:32,313,000-32,317,000" 5 0.15
run_region gnas_core   "chr20:58,899,000-58,904,000" 5 0.15
run_region h19_core    "chr11:1,999,000-2,003,000"   5 0.15
run_region snrpn_core  "chr15:24,954,000-24,957,000" 5 0.15

echo "" >> "$OUT"
echo "=== 全部 region RESULT 行 ===" >> "$OUT"
grep "RESULT" "$OUT" >> "$OUT"
echo "[done] $OUT"
