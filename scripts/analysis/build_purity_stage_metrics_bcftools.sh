#!/bin/bash
# build_purity_stage_metrics_bcftools.sh
#
# 目的：
# 1) 依 run_purity_and_standard_verification.sh 的 purity_status.tsv 彙整各 purity 三階段指標
# 2) Stage0 (ClairS) 使用 bcftools 計算 TP/FP/FN (限制於 truth BED)
# 3) 自動修正 pileup_filter.vcf 的 GQ header 型別 (Integer -> Float) 以避免 bcftools 解析失敗
#
# 使用方式：
#   bash scripts/analysis/build_purity_stage_metrics_bcftools.sh \
#     --sample HCC1395_DORADO \
#     --run-tag 20260213_dorado_purity_full

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../pipeline"
source "${PIPELINE_DIR}/config.sh"

SAMPLE="HCC1395"
RUN_TAG=""
OUTPUT_ROOT="/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup"
STATUS_TSV=""

usage() {
  cat <<EOF
Usage:
  $0 --sample NAME --run-tag TAG [--output-root DIR] [--status-tsv FILE]

Options:
  --sample NAME       樣本名稱（需存在於 pipeline config）
  --run-tag TAG       run_purity_and_standard_verification 的 run tag
  --output-root DIR   purity output root (預設: ${OUTPUT_ROOT})
  --status-tsv FILE   直接指定 status tsv（若提供則優先）
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2 ;;
    --run-tag) RUN_TAG="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --status-tsv) STATUS_TSV="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "${RUN_TAG}" ]]; then
  echo "[ERROR] --run-tag is required" >&2
  exit 1
fi

eval "$(get_sample_config "${SAMPLE}")"

if [[ -z "${STATUS_TSV}" ]]; then
  STATUS_TSV="${OUTPUT_ROOT}/${SAMPLE}/purity_runs/${RUN_TAG}/purity_status.tsv"
fi

if [[ ! -f "${STATUS_TSV}" ]]; then
  echo "[ERROR] status tsv not found: ${STATUS_TSV}" >&2
  exit 1
fi

RUN_DIR="$(dirname "${STATUS_TSV}")"
TMP_DIR="${RUN_DIR}/stage_metrics_tmp"
mkdir -p "${TMP_DIR}"

OUT_TSV="${RUN_DIR}/purity_stage_metrics.tsv"

echo "[INFO] SAMPLE=${SAMPLE}"
echo "[INFO] RUN_TAG=${RUN_TAG}"
echo "[INFO] STATUS_TSV=${STATUS_TSV}"
echo "[INFO] TRUTH_VCF=${TRUTH_VCF}"
echo "[INFO] TRUTH_BED=${TRUTH_BED}"

if [[ -z "${TRUTH_BED}" ]]; then
  echo "[ERROR] TRUTH_BED is empty for ${SAMPLE}. This script currently requires BED for stage0 evaluation." >&2
  exit 1
fi

TRUTH_KEYS="${TMP_DIR}/truth_keys.tsv"
"${BCFTOOLS}" view -R "${TRUTH_BED}" -v snps "${TRUTH_VCF}" | \
  "${BCFTOOLS}" query -f '%CHROM\t%POS\t%REF\t%ALT\n' | \
  sort -u > "${TRUTH_KEYS}"

TRUTH_TOTAL=$(wc -l < "${TRUTH_KEYS}")
echo "[INFO] TRUTH_TOTAL=${TRUTH_TOTAL}"

echo -e "subsample\tpurity\tmm_tags\tml_tags\tclairs_tp\tclairs_fp\tclairs_fn\tclairs_precision\tclairs_recall\tclairs_f1\tlongphase_tp\tlongphase_fp\tlongphase_fn\tlongphase_precision\tlongphase_recall\tlongphase_f1\ttp_regions_analyzed\tfp_regions_analyzed\tfinal_tp\tfinal_fp\tfinal_fn\tfinal_precision\tfinal_recall\tfinal_f1\tdelta_clairs_to_longphase\tdelta_longphase_to_final\tdelta_clairs_to_final\tnotes\tcleanup" > "${OUT_TSV}"

tail -n +2 "${STATUS_TSV}" | while IFS= read -r line; do
  run_tag=$(printf '%s\n' "${line}" | cut -f1)
  subsample=$(printf '%s\n' "${line}" | cut -f2)
  estimated_purity=$(printf '%s\n' "${line}" | cut -f3)
  tumor_bam=$(printf '%s\n' "${line}" | cut -f4)
  somatic_vcf=$(printf '%s\n' "${line}" | cut -f5)
  source_mm_tags=$(printf '%s\n' "${line}" | cut -f6)
  source_ml_tags=$(printf '%s\n' "${line}" | cut -f7)
  step01=$(printf '%s\n' "${line}" | cut -f8)
  step02=$(printf '%s\n' "${line}" | cut -f9)
  step03=$(printf '%s\n' "${line}" | cut -f10)
  cleanup=$(printf '%s\n' "${line}" | cut -f11)
  tp_regions=$(printf '%s\n' "${line}" | cut -f12)
  fp_regions=$(printf '%s\n' "${line}" | cut -f13)
  baseline_f1=$(printf '%s\n' "${line}" | cut -f14)
  filtered_f1=$(printf '%s\n' "${line}" | cut -f15)
  f1_delta=$(printf '%s\n' "${line}" | cut -f16)
  notes=$(printf '%s\n' "${line}" | cut -f17)
  output_dir=$(printf '%s\n' "${line}" | cut -f18)
  # sanitize ClairS vcf (GQ type patch) before bcftools parsing
  SANITIZED_VCF="${TMP_DIR}/${subsample}.sanitized.vcf.gz"
  if [[ ! -f "${SANITIZED_VCF}" ]]; then
    if [[ "${somatic_vcf}" == *.gz ]]; then
      CAT_CMD="zcat"
    else
      CAT_CMD="cat"
    fi
    "${CAT_CMD}" "${somatic_vcf}" | \
      sed 's/ID=GQ,Number=1,Type=Integer/ID=GQ,Number=1,Type=Float/' | \
      "${BCFTOOLS}" view -f PASS -Oz -o "${SANITIZED_VCF}"
    "${BCFTOOLS}" index -f "${SANITIZED_VCF}"
  fi

  CLAIRS_KEYS="${TMP_DIR}/${subsample}.clairs_keys.tsv"
  "${BCFTOOLS}" view -R "${TRUTH_BED}" -v snps "${SANITIZED_VCF}" | \
    "${BCFTOOLS}" query -f '%CHROM\t%POS\t%REF\t%ALT\n' | \
    sort -u > "${CLAIRS_KEYS}"

  clairs_tp=$(comm -12 "${CLAIRS_KEYS}" "${TRUTH_KEYS}" | wc -l)
  clairs_fp=$(comm -23 "${CLAIRS_KEYS}" "${TRUTH_KEYS}" | wc -l)
  clairs_fn=$((TRUTH_TOTAL - clairs_tp))

  clairs_precision=$(awk -v tp="${clairs_tp}" -v fp="${clairs_fp}" 'BEGIN{if(tp+fp==0) print 0; else printf "%.4f", tp/(tp+fp)}')
  clairs_recall=$(awk -v tp="${clairs_tp}" -v fn="${clairs_fn}" 'BEGIN{if(tp+fn==0) print 0; else printf "%.4f", tp/(tp+fn)}')
  clairs_f1=$(awk -v p="${clairs_precision}" -v r="${clairs_recall}" 'BEGIN{if(p+r==0) print 0; else printf "%.4f", 2*p*r/(p+r)}')

  # longphase stage
  if [[ -f "${output_dir}/longphase_s/variant_counts.txt" ]]; then
    # shellcheck source=/dev/null
    source "${output_dir}/longphase_s/variant_counts.txt"
    longphase_tp="${TP_COUNT}"
    longphase_fp="${FP_COUNT}"
    longphase_fn=$((TRUTH_TOTAL - longphase_tp))
  else
    longphase_tp=0
    longphase_fp=0
    longphase_fn="${TRUTH_TOTAL}"
  fi

  longphase_precision=$(awk -v tp="${longphase_tp}" -v fp="${longphase_fp}" 'BEGIN{if(tp+fp==0) print 0; else printf "%.4f", tp/(tp+fp)}')
  longphase_recall=$(awk -v tp="${longphase_tp}" -v fn="${longphase_fn}" 'BEGIN{if(tp+fn==0) print 0; else printf "%.4f", tp/(tp+fn)}')
  longphase_f1=$(awk -v p="${longphase_precision}" -v r="${longphase_recall}" 'BEGIN{if(p+r==0) print 0; else printf "%.4f", 2*p*r/(p+r)}')

  # final stage
  if [[ -f "${output_dir}/metrics.json" ]]; then
    read -r final_tp final_fp final_precision final_recall final_f1 <<< "$(python3 - "${output_dir}/metrics.json" <<'PY'
import json, sys
with open(sys.argv[1]) as f:
    m=json.load(f)
fl=m.get('filtered',{})
print(int(fl.get('tp',0)), int(fl.get('fp',0)), float(fl.get('precision',0)), float(fl.get('recall',0)), float(fl.get('f1',0)))
PY
)"
  else
    final_tp=0
    final_fp=0
    final_precision=0
    final_recall=0
    final_f1=0
  fi
  final_fn=$((TRUTH_TOTAL - final_tp))

  delta_clairs_to_longphase=$(awk -v a="${longphase_f1}" -v b="${clairs_f1}" 'BEGIN{printf "%.4f", a-b}')
  delta_longphase_to_final=$(awk -v a="${final_f1}" -v b="${longphase_f1}" 'BEGIN{printf "%.4f", a-b}')
  delta_clairs_to_final=$(awk -v a="${final_f1}" -v b="${clairs_f1}" 'BEGIN{printf "%.4f", a-b}')

  echo -e "${subsample}\t${estimated_purity}\t${source_mm_tags}\t${source_ml_tags}\t${clairs_tp}\t${clairs_fp}\t${clairs_fn}\t${clairs_precision}\t${clairs_recall}\t${clairs_f1}\t${longphase_tp}\t${longphase_fp}\t${longphase_fn}\t${longphase_precision}\t${longphase_recall}\t${longphase_f1}\t${tp_regions}\t${fp_regions}\t${final_tp}\t${final_fp}\t${final_fn}\t${final_precision}\t${final_recall}\t${final_f1}\t${delta_clairs_to_longphase}\t${delta_longphase_to_final}\t${delta_clairs_to_final}\t${notes}\t${cleanup}" >> "${OUT_TSV}"
done

echo "[INFO] Wrote: ${OUT_TSV}"
