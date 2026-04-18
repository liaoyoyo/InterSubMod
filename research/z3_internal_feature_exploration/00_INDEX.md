---
title: Z3 Internal Feature Exploration — Index
date: 2026-04-18
status: completed (NEGATIVE main + HCC1954-specific CONDITIONAL)
---

# Z3 內部特徵探索研究索引

## 動機

Zone-Aware Framework 確認 Z3（TO mode, LOH + extreme AF + HPFineNGroups≤1）跨樣本 TP rate σ=0.28（0.050–0.868），懷疑 Z3 內部仍有未探索的二階區分特徵。

## 執行步驟

| Step | 假說 | 腳本 | 結果 |
|------|------|------|------|
| 1 | H-Z3a/b/c：Z3 內有 AUC≥0.60 特徵 | `scripts/step1_z3_feature_auc_table.py` | NEGATIVE — 無特徵達 gate |
| 2.5 | H-Z3d：AF∈[0.4,0.6] germline FP 可用 CN × NG 分離 | `scripts/step2_5_af_germline_band.py` | NEGATIVE — 僅 1/7 樣本符合 |
| 3 | HCC1954 vs HCC1395 機制對比 | `scripts/step3_z3_mechanism_contrast.py` | POSITIVE 觀察 — HCC1954 amplicon-driven FP |

## 主要產出

- **報告**：`docs/experiments/in_progress/2026/04/20260418_Z3_Internal_Feature_Exploration_01.md`
- **數據**：`data/z3_feature_auc.tsv`, `data/z3_af_germline_band_*.tsv`, `data/z3_mechanism_*.tsv`
- **圖表**：
  - `figures/z3_feature_auc_heatmap.png`
  - `figures/z3_af_germline_band_heatmap.png`
  - `figures/z3_hcc1954_vs_hcc1395_mechanism.png`
  - `figures/z3_chr_distribution_contrast.png`

## 最終結論

1. **主結論 NEGATIVE**：Z3 內部無跨樣本可用的二階區分特徵；AF germline band 假說（H-Z3d）不成立
2. **HCC1954 特殊機制 POSITIVE（已知）**：Z3 FP 集中於 chr5/8/17（HER2/MYC amplicon），由 CNV artifact 驅動；已被 F pilot canonical filter 覆蓋
3. **對 5 研究目標**：僅 F1 目標 HCC1954-specific CONDITIONAL（無新增路徑）；其餘目標零衝擊
4. **Zone-Aware Framework 定位不變**：僅作 characterization annotation 使用
