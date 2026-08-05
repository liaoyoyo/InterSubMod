#!/bin/bash
# run_all_samples.sh — 對全 7 樣本 canonical paired_full longphase-s tagged BAM 跑
# 單次 ISM scan(tumor-only, distance matrix output) → compute_split_accounting.py → 切割總帳
#
# 流程(每樣本): ISM scan → driver(分群+Q1+分帶) → rmtree per-region matrices(保留 sig csv + records + result)
# 循序執行(每個 ISM 用 -j48 吃滿 48 核); disk-safe(逐樣本清理)。
set -uo pipefail

REPO=/big7_disk/liaoyoyo2001/InterSubMod
WT="$REPO/.claude/worktrees/ism-review-infra"
BIN="$WT/build/bin/inter_sub_mod"
REF=/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta
CANON=/big7_disk/liaoyoyo2001/big7_disk_output/canonical
SCRATCH=/big7_disk/liaoyoyo2001/ism_split_scan
RESULT="$REPO/docs/experiments/in_progress/2026/06/20260620_allsample_subcluster_split/results"
DRIVER="$REPO/scripts/analysis/subcluster_split_allsamples/compute_split_accounting.py"
export TMPDIR=/big7_disk/liaoyoyo2001/tmp
NPROC=46
THREADS=48
VTYPE="${VTYPE:-tp}"   # tp (default) 或 fp

mkdir -p "$RESULT" "$SCRATCH"

# 樣本清單 (小→大, 失敗影響最小先驗)
SAMPLES="HCC1937 HCC1954 HCC1395 HCC1395_DORADO COLO829 H1437 H2009"

echo "================ ALL-SAMPLE SUBCLUSTER SPLIT [VTYPE=$VTYPE] — $(date '+%F %T') ================"
echo "BIN=$BIN"
[[ -x "$BIN" ]] || { echo "FATAL: binary 不存在或不可執行: $BIN"; exit 1; }

for S in $SAMPLES; do
  echo ""
  echo "######## [$S] $(date '+%F %T') ########"
  PF=$(ls -d "$CANON/$S/paired_full/"*_complete_matrix 2>/dev/null | head -1)
  BAM="$PF/longphase_s/${S}_tagged.bam"
  TPVCF="$PF/longphase_s/filtered_snv_${VTYPE}.vcf.gz"
  OUT="$SCRATCH/${S}_${VTYPE}"
  # driver sample label: fp 加 _FP 後綴避免覆蓋 tp 結果
  if [[ "$VTYPE" == "fp" ]]; then SLAB="${S}_FP"; else SLAB="$S"; fi

  if [[ ! -f "$BAM" ]]; then echo "[$S] SKIP: tagged BAM 不存在 ($BAM)"; continue; fi
  if [[ ! -f "$TPVCF" ]]; then echo "[$S] SKIP: ${VTYPE} VCF 不存在 ($TPVCF)"; continue; fi

  NVAR=$(bcftools view -H "$TPVCF" 2>/dev/null | wc -l)
  echo "[$S] BAM=$BAM"
  echo "[$S] ${VTYPE^^} VCF variants=$NVAR"
  echo "[$S] 磁碟剩餘: $(df -h /big7_disk | tail -1 | awk '{print $4}')"

  # ---- Step 1: ISM scan (tumor-only, distance matrix output 預設 enabled) ----
  rm -rf "$OUT"; mkdir -p "$OUT"
  echo "[$S] [1/3] ISM scan (threads=$THREADS, w=5000, BERNOULLI, C_min=3, SKIP)..."
  T1=$(date +%s)
  "$BIN" -t "$BAM" -r "$REF" -v "$TPVCF" \
    -w 5000 -j "$THREADS" --distance-metric BERNOULLI --min-common-coverage 3 \
    --nan-distance-strategy SKIP --methyl-low 0.2 --methyl-high 0.8 \
    -o "$OUT" --log-level info > "$OUT/ism_run.log" 2>&1
  RC=$?
  T2=$(date +%s)
  if [[ $RC -ne 0 ]]; then echo "[$S] FATAL: ISM exit=$RC (見 $OUT/ism_run.log)"; tail -5 "$OUT/ism_run.log"; continue; fi
  NMAT=$(find "$OUT" -path "*distance/BERNOULLI/matrix.csv" 2>/dev/null | wc -l)
  echo "[$S] ISM done in $((T2-T1))s; matrices=$NMAT; sig csv=$(ls -la "$OUT/significance_summary.csv" 2>/dev/null | awk '{print $5}')"

  # ---- Step 2: driver (分群 + Q1 + 分帶) ----
  echo "[$S] [2/3] driver clustering + Q1 + bands (nproc=$NPROC)..."
  python3 "$DRIVER" "$SLAB" "$OUT" "$RESULT" "$NPROC"
  RC=$?
  if [[ $RC -ne 0 ]]; then echo "[$S] FATAL: driver exit=$RC"; continue; fi

  # ---- Step 3: 保留 sig csv + records + result; 清 per-region matrices ----
  echo "[$S] [3/3] 保留 significance_summary.csv, 清 per-region 目錄..."
  cp "$OUT/significance_summary.csv" "$RESULT/significance_summary_${SLAB}.csv" 2>/dev/null
  cp "$OUT/ism_run.log" "$RESULT/ism_run_${SLAB}.log" 2>/dev/null
  # 清掉龐大的 per-region 子目錄(chr*)，保留 scratch 根的 csv/log
  find "$OUT" -maxdepth 2 -type d -name "chr*" -exec rm -rf {} + 2>/dev/null
  rm -rf "$OUT"/*/chr* 2>/dev/null
  echo "[$S] DONE. result=$RESULT/result_${S}.json"
done

echo ""
echo "================ ALL DONE — $(date '+%F %T') ================"
echo "結果目錄: $RESULT"
ls -la "$RESULT"/result_*.json 2>/dev/null
