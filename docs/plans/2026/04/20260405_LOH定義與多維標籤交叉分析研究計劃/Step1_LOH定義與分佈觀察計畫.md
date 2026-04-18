<!--
建立時間: 2026-04-05 23:30
目標: Step 1 — LOH 雙定義交叉分析與 TP/FP 分佈觀察
處理範圍: ISM HP Imbalance vs LOH.bed 差異量化
關聯檔案:
  - docs/plans/2026/04/20260405_LOH定義與多維標籤交叉分析研究計劃/00_總覽與執行順序.md
  - docs/architecture/20260405_LOH定義與多維標籤交叉分析研究計劃/20260405_研究方向修正與觀察策略確認_01.md
-->

# Step 1: LOH 定義與分佈觀察計畫

## 背景

ISM 有兩種 LOH 相關的定義，從未被系統性比較：

| 指標 | 來源 | 判定方式 | 驗證狀態 |
|------|------|----------|----------|
| **HP Imbalance** (`Potential_LOH`) | ISM 內部計算 | HP_Ratio < 0.1 or > 0.9 | 已知受 self-phasing 污染 |
| **LOH.bed** (`to_loh_bed_hit`) | LongPhase-TO phasing 輸出 | 基因體區間 BED 格式 | SEQC2 Jaccard=0.8470 |

**術語**：本研究中 ISM 的 `Potential_LOH` 改稱「**HP Imbalance（HP 不平衡）**」。

## 資料來源

- **Master dataset**: `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz`
  - 748,391 rows × 116 columns
  - 關鍵欄位: `Potential_LOH`, `to_loh_bed_hit`, `HP_Ratio`, `HP1FamilyN`, `HP2FamilyN`, `truth_label`, `sample`, `mode`, `Coverage_Multiple`, `Coverage_Category`, `VerificationClass`, `NumReads`, `caller_af`
- **分析公用函式庫**: `scripts/analysis/observation_common.py`
  - `load_master_dataset()`, `setup_plot_style()`, `save_figure()`, `compute_auc()`, `compute_enrichment()`
  - `SAMPLE_ORDER`, `MODE_ORDER`, `OUTPUT_ROOT`, `TRUTH_PALETTE`, `MODE_PALETTE`
- **LongPhase 版本**: 所有分析以 LongPhase-TO **baseline** 為基準

## 子任務

### 1.1 LOH 雙定義交叉分析

**優先度**: P0（基礎）
**可並行**: 是（獨立於其他子任務）

**目標**: 量化 `Potential_LOH`（HP Imbalance）與 `to_loh_bed_hit`（LOH.bed）的一致性

**分析步驟**:

1. 載入 master dataset（TO mode only）
2. 建立 2×2 混淆矩陣：`Potential_LOH` × `to_loh_bed_hit`，per sample
3. 計算每個 sample 的：
   - Cohen's kappa（一致性）
   - Jaccard index
   - 敏感度（LOH.bed 為 reference 時的 HP Imbalance 敏感度）
   - 特異度
4. 繪製 HP_Ratio 連續分佈直方圖，按 `to_loh_bed_hit` 著色
5. 分析 HP_Ratio 閾值的 ROC（以 `to_loh_bed_hit` 為 ground truth）

**圖表規格**:

| 圖 | 內容 | 維度 |
|----|------|------|
| F1 | 2×2 confusion matrix heatmap | 7 samples panel (TO only) |
| F2 | HP_Ratio distribution, colored by LOH.bed | 7 samples × histogram |
| F3 | Cohen's kappa bar chart | 7 samples |
| F4 | HP_Ratio ROC curve | 7 samples overlay |

**輸出**:
- `{OUTPUT_ROOT}/20260406_loh_concordance/01_confusion_matrix.png`
- `{OUTPUT_ROOT}/20260406_loh_concordance/02_hp_ratio_distribution.png`
- `{OUTPUT_ROOT}/20260406_loh_concordance/03_kappa_barchart.png`
- `{OUTPUT_ROOT}/20260406_loh_concordance/04_hp_ratio_roc.png`
- `{OUTPUT_ROOT}/20260406_loh_concordance/concordance_statistics.tsv`

### 1.2 四象限 TP/FP 分佈

**優先度**: P0（基礎）
**可並行**: 是

**目標**: 在每個 LOH 象限中分析 TP/FP 的分佈與比例

**四象限定義**:
```python
df['loh_quadrant'] = 'Q4_neither'  # default
mask_ism = df['Potential_LOH'] == True
mask_bed = df['to_loh_bed_hit'] == True
df.loc[mask_ism & mask_bed, 'loh_quadrant'] = 'Q1_both'
df.loc[mask_ism & ~mask_bed, 'loh_quadrant'] = 'Q2_ism_only'
df.loc[~mask_ism & mask_bed, 'loh_quadrant'] = 'Q3_bed_only'
```

**分析步驟**:
1. 計算每象限 × 每 sample 的 TP count, FP count, FP rate
2. 計算 FP enrichment ratio（vs 全域 FP rate）
3. Fisher exact test（每象限 TP/FP vs 其餘）

**圖表規格**:

| 圖 | 內容 | 維度 |
|----|------|------|
| F5 | Stacked bar: TP/FP per quadrant | 7 samples panel |
| F6 | FP rate heatmap (sample × quadrant) | 7×4 heatmap |
| F7 | FP enrichment ratio bar chart | 7 samples × 4 quadrants |

**輸出**:
- `{OUTPUT_ROOT}/20260406_loh_concordance/05_quadrant_tp_fp_stacked.png`
- `{OUTPUT_ROOT}/20260406_loh_concordance/06_fp_rate_heatmap.png`
- `{OUTPUT_ROOT}/20260406_loh_concordance/07_fp_enrichment.png`
- `{OUTPUT_ROOT}/20260406_loh_concordance/quadrant_statistics.tsv`

### 1.3 差異區域生物學解釋

**優先度**: P1（依賴 1.2）
**可並行**: 否（需 1.2 結果）

**目標**: Q2（ISM-only HP Imbalance）和 Q3（LOH.bed-only）的特徵差異與生物學解釋

**Q2 假說**（ISM-only，HP Imbalance 但不在 LOH.bed 內）:
- Self-phasing artifact（ALT reads 偏向一個 HP）
- 低 read count 的隨機波動（NumReads < 20）
- 真正的 allelic imbalance（variant 確實在一個 haplotype，但不在基因體 LOH 區域）

**Q3 假說**（LOH.bed-only，在 LOH.bed 但 HP 平衡）:
- LOH.bed 覆蓋大段 region，內部有少數 HP 平衡位點
- Copy-neutral LOH（coverage 正常但失去 heterozygosity）
- Phase block boundary

**分析步驟**:
1. 比較 Q2 vs Q3 的特徵分佈：
   - `NumReads`（Q2 預期較低 → 隨機波動）
   - `Coverage_Multiple`（Q3+正常 coverage → cnLOH 候選）
   - `caller_af`（Q2 高 AF 可能是 self-phasing）
   - `VerificationClass` 分佈
   - `HP_Ratio` 連續分佈（Q2 接近 0/1，Q3 接近 0.5）
2. Q2 中 `NumReads < 20` 的比例 → 量化隨機波動貢獻
3. Q3 中 `Coverage_Multiple ∈ [0.8, 1.2]` 的比例 → 量化 cnLOH 候選

**圖表規格**:

| 圖 | 內容 |
|----|------|
| F8 | Violin plots: NumReads by quadrant (TO only) |
| F9 | Violin plots: Coverage_Multiple by quadrant |
| F10 | Violin plots: caller_af by quadrant |
| F11 | VerificationClass stacked bar by quadrant |

### 1.4 SEQC2 驗證（HCC1395 only）

**優先度**: P1
**可並行**: 是

**目標**: 用 SEQC2 truth set 驗證兩種 LOH 定義的準確性

**已有基礎**: O15 報告 LOH.bed Jaccard=0.8470

**分析步驟**:
1. 取 HCC1395 TO mode 數據
2. 以 SEQC2 truth label 為 reference
3. 計算 HP_Ratio 作為 LOH predictor 的 AUC
4. 掃描 HP_Ratio 閾值 (0.05-0.20 and 0.80-0.95)，找最佳 F1
5. 與 LOH.bed 的 Jaccard 比較

**輸出**: `{OUTPUT_ROOT}/20260406_loh_concordance/seqc2_validation.tsv`

## 腳本模板

```python
#!/usr/bin/env python3
"""LOH Definition Concordance Analysis: HP Imbalance vs LOH.bed"""

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from observation_common import *

OUT_DIR = ensure_dir(OUTPUT_ROOT / "20260406_loh_concordance")

def main():
    df = load_master_dataset()
    df_to = df[df['mode'] == 'to'].copy()
    
    # 1.1 Concordance
    # ...
    
    # 1.2 Quadrants
    # ...
    
    # 1.3 Feature comparison (if 1.2 done)
    # ...
    
    # 1.4 SEQC2 (HCC1395 only)
    # ...

if __name__ == "__main__":
    main()
```

## 驗證清單

- [ ] HP Imbalance rate（TO）與 O3 報告一致（~60%）
- [ ] LOH.bed Jaccard（HCC1395）重現 O15 的 0.8470
- [ ] 四象限 TP+FP 合計 = 每 sample 總數
- [ ] 所有圖表跨 7 samples 完成
- [ ] 不一致區域（Q2, Q3）的生物學解釋有數據支持

## 後續銜接

Step 1 完成後：
1. 向使用者報告結果，確認是否啟動 Step 2（HP × ALT/REF 提取）
2. Step 3.1（PERMANOVA 審計）可與 Step 1 並行
3. Q3（LOH.bed-only + 正常 coverage）的發現直接餵入 Step 3.2（cnLOH 分析）
