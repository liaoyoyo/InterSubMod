#!/bin/bash
# T1.2-F1 全基因組 priority bug audit 一鍵腳本
# 前置：vote_dump_{baseline,v3f,v5}_genome.tsv.gz 已在當前目錄
# 使用：cd InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit && bash run_t12_f1_genome_analysis.sh

set -euo pipefail

REPO=/big7_disk/liaoyoyo2001/InterSubMod
DUMP_DIR=${REPO}/research/v5_provenance_followup/T1_2_read_level_audit
OUT_MD=${DUMP_DIR}/T1_2_F1_genome_wide_audit.md

echo "[T1.2-F1] === 全基因組 priority bug audit 啟動 $(date) ==="
echo "[1/4] 確認輸入 dump 檔案"
ls -lh ${DUMP_DIR}/vote_dump_baseline_genome.tsv.gz \
       ${DUMP_DIR}/vote_dump_v3f_genome.tsv.gz \
       ${DUMP_DIR}/vote_dump_v5_genome.tsv.gz

echo "[2/4] 跑 4 路徑 + per-chr + Layer 1.5 trigger detection"
cd ${REPO}
python3 scripts/analysis/v5_provenance_vote_audit.py \
    --dump-dir ${DUMP_DIR} \
    --suffix genome \
    --out ${OUT_MD}

echo "[3/4] 報告產出在 ${OUT_MD}"
ls -lh ${OUT_MD}

echo "[4/4] 確認 .gitignore 排除 dump"
cat ${DUMP_DIR}/.gitignore | grep -E '^vote_dump'

echo "[T1.2-F1] === 完成 $(date) ==="
echo
echo "下一步建議："
echo "  cat ${OUT_MD}                                         # 看完整報告"
echo "  git add InterSubMod/scripts/analysis/v5_provenance_vote_audit.py"
echo "  git add InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/T1_2_F1_genome_wide_audit.md"
echo "  git add InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/run_t12_f1_genome_analysis.sh"
echo "  git status -s | grep vote_dump_  # 應該完全沒輸出（gitignore 生效）"
