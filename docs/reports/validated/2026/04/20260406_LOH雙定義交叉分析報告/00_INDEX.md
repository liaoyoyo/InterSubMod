<!--
建立時間: 2026-04-06 02:30
目標: LOH 雙定義交叉分析 Wave 1+2 完整報告索引
處理範圍: LOH 雙定義一致性、四象限 TP/FP、PERMANOVA 適用性、Q2 解釋、SEQC2 驗證、Allele AUC
關聯檔案:
  - docs/plans/2026/04/20260405_LOH定義與多維標籤交叉分析研究計劃/00_總覽與執行順序.md
  - docs/architecture/20260405_LOH定義與多維標籤交叉分析研究計劃/20260405_研究方向修正與觀察策略確認_01.md
-->

# LOH 雙定義交叉分析報告

**日期**: 2026-04-06 | **LongPhase**: Baseline | **數據**: Master dataset 748,391 rows × 116 cols

## 報告結構

| # | 章節 | 對應 Step | 圖表 | 統計表 |
|---|------|----------|------|--------|
| 01 | [LOH 雙定義一致性與四象限](01_LOH雙定義一致性與四象限.md) | 1.1 + 1.2 | 7 張 | 3 個 |
| 02 | [PERMANOVA 適用性審計](02_PERMANOVA適用性審計.md) | 3.1 | 8 張 | 2 個 |
| 03 | [Q2 差異區域解釋](03_Q2差異區域解釋.md) | 1.3 | 4 張 | 2 個 |
| 04 | [SEQC2 驗證與 F1 影響](04_SEQC2驗證與F1影響.md) | 1.4 | 5 張 | 4 個 |
| 05 | [Allele AUC 與 F1 影響](05_Allele_AUC與F1影響.md) | 3.3 | 6 張 | 5 個 |
| 06 | [綜合判定與後續方向](06_綜合判定與後續方向.md) | — | — | — |
| 07 | [Wave 3: 非 LOH 區分力與特徵方向](07_Wave3_非LOH區分力與特徵方向.md) | M2+M3 | 11 張 | 6 個 |
| 08 | [Wave 3: AlleleDelta 深入與多特徵組合](08_Wave3_AlleleDelta深入與多特徵組合.md) | M4 | 5 張 | 5 個 |
| 09 | [Wave 3: Bug Fix + 肉眼檢視 + ISM 改進](09_Wave3_肉眼檢視文件與ISM改進.md) | M1+M5+M6 | 120 張 | 1 個 |
| **10** | [**完整圖表比較檢視**](10_完整圖表比較檢視.md) | 全 Wave | **166 張全內嵌** | 40 對照 |

## 核心判定速查

| # | 判定 | 置信度 |
|---|------|--------|
| J1 | HP Imbalance 是 LOH.bed 的超集 (Sensitivity 99.7-100%) | 極高 |
| J2 | HP PERMANOVA 在 LOH 內不可用 (Valid rate 5-6%) | 確定 |
| J6 | Q2 是非 LOH 區域的 HP phasing 偏差 (AF=0.47, d=-1.04) | 高 |
| J7 | LOH.bed 與 SEQC2 一致性極高 (Jaccard=0.928) | 極高 |
| J9 | **LOH 區域不可作為 variant filter** (10/10 情境 FAIL TP loss ≤ 2%) | **確定** |
| J10 | AlleleDelta 在 LOH 內有真實但微弱的區分力 (AUC=0.556, confound-free) | 中 |
| J11 | Non-LOH 區分力同樣有限 (max AUC=0.643, read count proxy) | 高 |
| J13 | **LOH 多特徵組合不可行** (Voting AUC=0.577 < 0.58) | **確定** |
| J14 | cnLOH PairwiseMeanDist 0.587 是 Simpson's Paradox (mean per-sample=0.50) | 確定 |
| J15 | CramersV 在 LOH 被 NumReads confound (0.511→0.464) | 高 |
| J16 | AlleleDelta 是 LOH 唯一 confound-free 信號 (但 AUC=0.556 不足) | 確定 |

## 戰略結論

> **LOH 的價值在 stratification（分層分析），而非 filtering（篩選移除）。**

## 數據來源

| 資料 | 路徑 |
|------|------|
| Master dataset | `big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz` |
| SEQC2 LOH benchmark | `/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed` |
| 分析腳本 (W1+W2) | `scripts/analysis/build_loh_concordance_analysis.py`, `build_permanova_loh_audit.py`, `build_loh_quadrant_explanation.py`, `build_seqc2_loh_validation.py`, `build_allele_loh_auc_f1_analysis.py` |
| 分析腳本 (W3) | `build_non_loh_discrimination.py`, `build_feature_direction_map.py`, `build_allele_deep_dive.py`, `build_visual_inspection_report.py` |
| 圖表原始輸出 | `big7_disk_output/synthesis/observation_workspaces/20260406_*`, `20260408_*` |
