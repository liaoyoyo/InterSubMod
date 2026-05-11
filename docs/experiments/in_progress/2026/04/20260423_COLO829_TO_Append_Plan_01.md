---
title: COLO829 TO append + Phase 3 synthesis rerun — 腳本使用說明
date: 2026-04-23
status: in_progress
owner: InterSubMod Research
related:
  - scripts/analysis/20260423_colo829_to_append_and_phase3_rerun.sh
  - scripts/analysis/20260423_phase3_synthesis.py
  - docs/experiments/in_progress/2026/04/20260423_Phase3_Cross_Sample_S1S7_Synthesis_01.md
  - research/ng_kde_rescaling/data/merged_7samples_paired_full_plus_7to_full_20260423.tsv.gz
---

# COLO829 TO append + Phase 3 synthesis rerun

## 用途

COLO829 TO full pipeline（`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260423_colo829_to_pilot/`）從 step01 ClairS-TO 起獨立重跑；pipeline 完成後使用本腳本：

1. 驗證 step05 intersubmod_tp/fp `significance_summary.csv` 存在且非空
2. 驗證新 KDE binary 輸出（`Diploid_Coverage_Used` 欄存在、非 NaN）
3. 讀取 COLO829 TO TP+FP，加 sample/mode/sample_order 等欄位
4. 與現有 `merged_7samples_paired_full_plus_hcc1395_to.tsv.gz` 做 row-append（欄位以 union 對齊，缺欄位填 NaN）
5. 輸出新 master `merged_7samples_paired_full_plus_7to_full_20260423.tsv.gz`
6. 執行 `20260423_phase3_synthesis.py`（臨時以 symlink 方式 swap-in 新 master，跑完自動還原原檔）
7. 列印 COLO829 TO 的 S1-S7 結果 + 與其他 6 TO 對照

## 呼叫方式

```bash
# 僅檢查 pipeline 輸出是否到位（不做任何修改）
bash scripts/analysis/20260423_colo829_to_append_and_phase3_rerun.sh --check-only

# 完整執行（append + phase3 rerun）
bash scripts/analysis/20260423_colo829_to_append_and_phase3_rerun.sh
```

腳本 **idempotent**：
- 若新 master 已存在會直接覆蓋
- 若 COLO829/to_pileup 已在 base master 會先移除再 append（避免重複）
- Step 4 phase3 swap 以 `trap` 保證原檔還原（成功/失敗皆還原）
- 會產生時間戳 `.bak.colo829_append_YYYYMMDD_HHMMSS` 備份

## 預期輸出

| 檔案 | 說明 |
|------|------|
| `research/ng_kde_rescaling/data/merged_7samples_paired_full_plus_7to_full_20260423.tsv.gz` | 新 master (~18 MB)，748,676 rows，14 個 (sample, mode) 組合 |
| `docs/experiments/in_progress/2026/04/figures/20260423_phase3_synthesis/summary.json` | Phase 3 最終 verdict（S3/S5 等） |
| `.../s1s7_per_sample.tsv` | 112 rows（7 TO + 7 paired × 8 schemes） |
| `.../s1s7_heatmap_tp_rate.png` | 熱圖含 n 標註 |
| `.../s1s7_heatmap_fold.png` | Fold vs baseline 熱圖 |
| `.../s3_cross_sample_wilcoxon.json` | S3 Wilcoxon p-value |
| `.../ng_gap_aggregate.json` | Thread D NG=2 gap（B1/B2 整合） |

## 初次執行結果（2026-04-23 09:30）

COLO829 TO pipeline 提前於 05:31 完成（非原預估 4-6 hr），腳本隨即端到端執行成功：

- TP rows=33,089; FP rows=17,528；KDE baseline=29.0×（與 KDE_Fix_Acceptance 表 §5.0 COLO829 29× 一致）
- merged master: 748,676 rows (base 698,059 + COLO829 TO 50,617)，14 個 (sample, mode) 齊全
- Phase 3 S3 Wilcoxon（TO, greater）p=**0.0625**（升格，6 樣本時為 p=0.125）
- S3 Verdict: INCONCLUSIVE → **CONDITIONAL_POSITIVE**
- S5 Verdict: **POSITIVE** 維持（median 0.876, 7/7 涵蓋）
- COLO829 TO 關鍵數據：baseline TP=0.654，S1 TP=0.656（與 baseline 差異微小，與 B5 結論「S1 fold 作 cross-sample metric 需樣本量閾值」一致）、S3 TP=0.718 (n=39, sparse)、S6/S7 n=10/20 過小

## 後續動作

- 將 Phase 3 synthesis doc 的 "6/6 TO" 宣告升級為 "7/7 TO"，並加入 S3 Wilcoxon p=0.0625 記錄
- COLO829 S1 TP rate 0.656 與其他 TO 0.876 median 差距大，呼應 B5 S1 fold=0.59 反轉分析；可納入 C2 archive TO 其他 6 樣本 rerun 後的整合視圖
