<!--
建立時間: 2026-04-06
目標: Wave 3 M1+M5+M6 結果 — Bug fix、ISM 改進方案、肉眼檢視文件
處理範圍: w2a_04 bug 修正、ISM 程式碼改進方案調整、40 位點代表性肉眼檢視
關聯檔案:
  - scripts/analysis/build_allele_loh_auc_f1_analysis.py (M1)
  - scripts/analysis/build_visual_inspection_report.py (M6)
  - big7_disk_output/synthesis/observation_workspaces/20260406_visual_inspection/
-->

# Wave 3: Bug Fix + 肉眼檢視文件 + ISM 改進方案

## M1: Bug Fix — w2a_04 圖表修正

### 問題

`build_allele_loh_auc_f1_analysis.py` 的 `fig04_valid_rate_comparison()` 使用 `np.isfinite(score)` 計算 Delta 欄位的 valid rate — 但 Delta 欄位 always present (非 NaN)，導致顯示 100% valid rate。

### 修正

改用 `LabelHPPermanovaValid` / `LabelAllelePermanovaValid` 布林欄位計算。修正後：

| 指標 | 修正前 | 修正後 | Wave 1 參考 |
|------|--------|--------|------------|
| TO LOH HP valid | 100% | **5.4%** | 5.4% ✓ |
| TO LOH Allele valid | 100% | **58.6%** | 58.6% ✓ |

### 修正位置

- `build_allele_loh_auc_f1_analysis.py` ~line 113-124 (整體 AUC matrix)
- `build_allele_loh_auc_f1_analysis.py` ~line 157-168 (per-sample section)

---

## M5: ISM 程式碼改進方案（修訂版）

### 前提條件變更

原計劃 M5 依賴 M4 結果。M4 判定：

- 多特徵組合 AUC = 0.577 < 0.58 → **NOT FEASIBLE**
- cnLOH PairwiseMeanDist per-sample 驗證 **FAIL**
- AlleleDelta AUC = 0.556 — 真實但太弱

### Code-1: LOH-Aware Dimension Selection — 降級為 P2

**原優先度 P0 → 降為 P2**

理由：M4 確認 Allele 維度在 LOH 內最高 AUC 僅 0.556，即使切換到 Allele 維度，區分力仍不足。切換邏輯的工程成本不值得這麼微弱的收益。

**保留為觀察性質**：在 diagnostic output 中記錄 LOH-aware 資訊（已有 `is_potential_loh()` flag），但不影響 scoring。

### Code-2: QualityScore TO 模式重設計 — 維持 P1

修正方向調整：

| 原方案 | 修正後方案 |
|--------|-----------|
| LOH 區域使用 Allele-based components | LOH 區域**取消所有反向特徵的貢獻**（caller_af、HPMergedDelta） |
| 加入 AlleleDelta 作為 LOH QS 組件 | 僅作微調（+0.004 AUC 不值得新增組件） |

核心修正：**消除反向特徵的負面影響**，而非加入新的正向特徵。

### Code-3: PERMANOVA Valid Rate 診斷輸出 — 維持 P2

不變，仍建議在 batch 結束時輸出 valid rate 統計。

### 不建議的修改（更新）

| 不建議 | 理由 | 新增/原有 |
|--------|------|----------|
| 降低 `min_reads_per_group` | 統計學底線 | 原有 |
| LOH-based variant filtering | J9 確定性關閉 | 原有 |
| LOH penalty 重新啟用 | O1-O10 已確認反向 | 原有 |
| **LOH-Aware dimension switching** | **M4: Allele AUC=0.556 不足** | **新增** |
| **多特徵組合 scoring** | **M4: 最佳 combo AUC=0.577 < 0.58** | **新增** |

---

## M6: 肉眼檢視文件

### 概述

40 個代表性位點（HCC1395 TO mode），按 4 個類別 × TP/FP × 5 抽樣。每個位點生成 3 張圖 + 文字摘要，共 **120 張圖 + 40 份摘要 + 4 份比較模板**。

### 抽樣分類

| 類別 | 條件 | TP | FP | 目的 |
|------|------|----|----|------|
| LOH + AlleleDelta 高 | Q1 + AlleleDelta > 75th | 5 | 5 | LOH 區域最有區分力的位點 |
| LOH + AlleleDelta 低 | Q1 + AlleleDelta < 25th | 5 | 5 | LOH 區域 hopeless 位點 |
| Non-LOH + HP 顯著 | Q4 + HP sig = True | 5 | 5 | Non-LOH 區域 HP 有效的位點 |
| Non-LOH + 無信號 | Q4 + HP/Allele sig = False | 5 | 5 | 最難區分的 non-LOH 位點 |

### 排除條件

- NumReads 限制 [30, 100]
- NumCpGs ≥ 10
- **排除著絲點與端粒區域**（hg38 centromere/telomere 座標）

### 每位點產出

| 圖表 | 說明 |
|------|------|
| `v1_methylation_heatmap.png` | Read × CpG 甲基化矩陣，按 HP tag 排序分色 |
| `v2_distance_heatmap.png` | Bernoulli 距離矩陣熱圖 + 階層聚類 dendrogram |
| `v3_read_structure.png` | HP tag × ALT/REF 分佈 + strand 結構 |
| `summary.md` | 位點基本資訊 + PERMANOVA 統計 + 特徵值 |

### IGV 截圖

IGV batch script 已生成 (`igv_batch.txt`)，包含 40 個位點的自動截圖指令。使用者可手動執行：

```bash
# 使用現有 X11 display
/bip7_disk/IGV_Linux_2.11.1/igv.sh --batch /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260406_visual_inspection/igv_batch.txt
```

### 比較模板

每個分類目錄下的 `comparison.md` 提供結構化的肉眼觀察紀錄模板：
- 甲基化模式差異
- HP / Allele 分佈差異
- IGV 讀段特徵
- 總結判斷

### 輸出路徑

```
big7_disk_output/synthesis/observation_workspaces/20260406_visual_inspection/
├── 00_inspection_index.md
├── sampled_positions.tsv
├── igv_batch.txt
├── igv_locus_list.tsv
├── igv_screenshots/
├── loh_high_allele/     (10 sites: 5 TP + 5 FP)
├── loh_low_allele/      (10 sites: 5 TP + 5 FP)
├── nonloh_hp_sig/       (10 sites: 5 TP + 5 FP)
└── nonloh_no_signal/    (10 sites: 5 TP + 5 FP)
```
