<!--
建立時間: 2026-03-27 00:00
目標: LOH Round 2 執行規格，聚焦在「有效 LOH evidence 的品質條件」
處理範圍: effective_hp support 分層分析、HP0 來源分析
關聯檔案:
  - docs/reports/validated/2026/03/20260327_LOH_round1_cross_sample_audit_01.md
  - docs/architecture/20260327_InterSubMod研究願景定錨_01.md
  - docs/plans/2026/03/20260326_LOH盤點執行規格_01.md
-->

# LOH Round 2 執行規格

> **接續**：Round 1 完成了「LOH-like 分佈與 TP/FP 差異」的診斷底圖。
> **Round 2 核心問題**：在什麼品質條件下，LOH evidence 才是可信的 evidence panel 成員？

---

## 一、Round 1 確認的前提

| 確認事項 | 結論 |
|---------|------|
| `HP_Ratio = 0.5` 的語意 | 不代表「平衡」，可能是 `effective_hp_reads = 0` |
| `HP1FamilyN` 的定義 | HP1-family read 數（HP1 + HP1-1），C++ LabelTest Stage 1 計算 |
| `HP2FamilyN` 的定義 | HP2-family read 數（HP2 + HP2-1），同上 |
| `NHP0` 的定義 | HP0（unphased）read 數，在 LOH ratio 計算中被排除 |
| `NHP3` 的定義 | HP3 read 數，同樣被排除 |
| LOH-like 在 paired | 整體 FP/TP enrichment 1.194×，但 sample heterogeneity 大 |
| LOH-like 在 TO | ~~TP/FP 幾乎同等 LOH-like（0.912×）~~ → **TP 富集 0.805×**（[修正 2026-03-30] HP integer tag bug），TO LOH 是 TP marker 而非 FP marker |
| TO effective_hp_reads | 大量 region `< 10`（COLO829 TP: 67%），LOH 訊號建立在弱 HP 支持上 |

**Round 2 的基本假說**：排除 `effective_hp_reads < 10` 的雜訊後，`LOH-like` 在 paired 端是否能展現更穩定、可泛化的 FP enrichment？

---

## 二、Core 1 — effective_hp support 分層分析

### 目標

確認 FP/TP enrichment 是否在 high-support tier 顯著且穩定提升。

### 分層定義

| Tier | 條件 | 預期狀況 |
|------|------|---------|
| Tier A | `effective_hp_reads ≥ 30` | 可信 LOH |
| Tier B | `effective_hp_reads 10–29` | 謹慎解讀 |
| Tier C | `effective_hp_reads < 10` | 幾乎不可信（Round 1 已確認為雜訊） |

### 分析內容

1. **paired 端**：每個 tier 內重算各樣本的 FP/TP LOH-like enrichment
   - 確認 Tier A 的 enrichment 是否穩定（7 樣本都 > 1.0×）
   - 確認 Tier A 的 enrichment 是否比 overall 更高（現在整體只有 1.194×）
   - 特別關注 HCC1954、H2009、H1437（Round 1 的高 enrichment 樣本）

2. **TO 端**：同樣分層，確認 Tier A 的 TO 中 LOH-like 是否仍然 TP/FP 同等高
   - 若 TO Tier A 也沒有 discriminability，進一步確認 TO 需要不同的分析框架

3. **paired vs TO 的 Tier A 比較**：同樣本同位點，Tier A LOH-like 的一致性

### 輸出格式

```
loh_tier_enrichment_by_sample_mode.tsv
  columns: sample, mode, tier, tp_total, fp_total, tp_loh_like_n, fp_loh_like_n,
           tp_loh_like_frac, fp_loh_like_frac, fp_tp_enrichment
```

---

## 三、Core 2 — HP0 來源分析

### 目標

確認 TO 中大量 HP0 的來源，與判斷 HP0 在 LOH-like region 的意義。

### 背景

C++ LabelTest 實作了 **Stage 4: Unassigned Affinity Test**：

- HP0 和 HP3 reads 被收集（excluded from HP ratio）
- 計算每個 unassigned read 到 HP1-family 和 HP2-family 的平均距離
- 給出 `affinity_score`：`(d_hp1 - d_hp2) / (d_hp1 + d_hp2)`
  - 正值表示 HP0 偏向 HP2-family
  - 負值表示 HP0 偏向 HP1-family

目前 Round 1 只用了 `NHP0` 和 `NHP3` 的數量（當 QC 欄位），未分析 affinity score 的分佈意義。

### 分析內容

**分析 1：LOH-like vs non-LOH-like 的 HP0 ratio 比較**

- 在 paired 和 TO 中分別比較 LOH-like region 的 `hp0_ratio` 分佈
- 假說：若 LOH 是真實的 one-sided phasing，HP0 應該不會特別高；若 LOH 是 phasing 失敗的假象，HP0 應該偏高

**分析 2：同位點 paired vs TO 的 HP0 ratio 差異**

- 利用 `same_locus_compare.tsv` 框架
- 確認 TO 的 HP0 是否系統性高於 paired 同位點
- 預期結論：TO phasing 在 LOH 區域更脆弱（longphase-to 缺少 paired 的 phasing constraints）

**分析 3：Stage 4 affinity score 分佈分析**

- 讀取 output summary 中的 `UnassignedAffinityScore` 欄位（如已輸出）
  - 若尚未在 output 中，需確認是否已在 summary TSV 中存在
- 分析：HP0 reads 的 affinity score 分佈在 LOH-like vs non-LOH-like region 是否不同
- 意義：若 HP0 在 LOH-like region 有強 affinity（傾向某一 family），可能補充 LOH evidence

**分析 4：TO 中高 HP0 的成因分類**

建立初步分類假說並嘗試驗證：

| 假說 | 驗證方式 |
|-----|---------|
| phasing 失敗（reads 無法被 longphase-to 分配） | 高 HP0 region 的 PS availability 低 |
| 真實 unphased reads（例如在 LOH 區域 phaser 無法決定） | 高 HP0 region 的 `core_loh_like` 也高，且同位點 paired 的 HP0 低 |
| 低 coverage 所致 | 高 HP0 region 的 `NumReads` 也低 |

### 輸出格式

```
hp0_loh_association_by_sample_mode.tsv
  columns: sample, mode, truth, loh_like_group, hp0_ratio_median, hp0_ratio_q25, hp0_ratio_q75,
           affinity_score_mean, affinity_score_std, n_regions

hp0_same_locus_delta.tsv
  columns: sample, locus, paired_hp0_ratio, to_hp0_ratio, hp0_delta,
           paired_loh_like, to_loh_like, concordance
```

---

## 四、附帶低優先項目（不在 Round 2 主範圍）

### longphase-s vs longphase-to LOH 一致性

- 同樣本同位點，兩種 phasing 下的 `LOH_Subtype` 是否一致
- **優先度**：低，等 Core 1/2 完成後再評估
- 說明：Round 1 已用 `paired` 對應 longphase-s、`TO` 對應 longphase-to，同位點比較已有框架

### LOH × methylation 聯合分析（暫緩至 Round 3）

- 目標：LOH-like region 的甲基分佈是否與 non-LOH-like 有差異？
- **暫緩原因**：需要在 support quality 分層確立後（Round 2），才能有意義地分析甲基差異
- Round 2 規格設計時預留 methylation 欄位的 join 介面（`methylation.csv` → `all_region_rows.tsv`）

---

## 五、資料欄位 Data Dictionary（補充說明）

此節補充 Round 1 中定義不夠清楚的欄位：

| 欄位 | 完整定義 | 計算來源 |
|------|---------|---------|
| `HP1FamilyN` | HP1-family read 數 = HP1 + HP1-1 reads | C++ `LabelTest::hp_to_merged_labels()` Stage 1，`n_hp1_family` |
| `HP2FamilyN` | HP2-family read 數 = HP2 + HP2-1 reads | 同上，`n_hp2_family` |
| `NHP0` | HP0（unphased）read 數，未被納入 HP ratio 計算 | C++ `n_hp0`（或等效欄位） |
| `NHP3` | HP3 read 數，未被納入 HP ratio 計算 | C++ `n_hp3` |
| `HP_Ratio` | 原始輸出（含 Laplace smoothing），**不能直接使用** | C++ `compute_hp_ratio(hp1_family_n, hp2_family_n)` |
| `effective_hp_reads` | Round 2 計算：`HP1FamilyN + HP2FamilyN`，排除 HP0/HP3 | Python 重算 |
| `hp_ratio_core` | Round 2 主定義：`HP1FamilyN / effective_hp_reads`，`NaN` when `effective_hp_reads = 0` | Python 重算 |
| `core_loh_like` | `effective_hp_reads > 0 AND (hp_ratio_core < 0.1 OR hp_ratio_core > 0.9)` | Round 1 定義，Round 2 繼承 |

---

## 六、執行前置確認

1. **確認 Stage 4 affinity score 是否已在 summary TSV 輸出**
   - 若無此欄位，Round 2 分析 3 需要從 per-region output 或重跑程式取得
2. **確認 `same_locus_compare.tsv` 的欄位是否包含 HP0 相關欄位**
   - 若無，需在 Round 2 腳本中從 `all_region_rows.tsv.gz` join
3. **確認 Round 2 的輸出 workspace 路徑**
   - 建議：`output/synthesis/observation_workspaces/20260327_loh_round2_support_stratification/`

---

## 七、Round 2 與研究願景的對應

| Round 2 分析 | 對應研究願景 |
|-------------|------------|
| effective_hp support 分層 | 目標 1 的可信 evidence 條件定義；目標 5 的 LOH feature 可用性確認 |
| HP0 來源分析 | 目標 4 的 TO phasing 品質理解；evidence panel 的 QC 機制 |
| Data dictionary 補充 | 基礎文件完整性，確保後續分析不會誤解欄位定義 |

---

## 八、Round 2 不應過度宣稱的事情

本輪**不能**直接宣稱：
- `LOH + support tier` 已可提升 TO 的 F1
- HP0 分析結果可直接當過濾條件

本輪**可以**穩定宣稱（完成後）：
- 哪個 support tier 的 LOH evidence 在 paired 端是穩定可信的
- TO 中高 HP0 的主要來源假說（phasing 失敗 vs 真實 unphased）
- Stage 4 affinity score 是否有補充 LOH evidence 的潛力
