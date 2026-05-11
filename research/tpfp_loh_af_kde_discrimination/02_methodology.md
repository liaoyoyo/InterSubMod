---
title: 方法學 — 五維切分、Scheme 定義、統計指標
status: in_progress
last_updated: 2026-04-22
---

# 02. 方法學

## 2.1 五個切分維度

這五個維度選擇 **以生物假設為先、資料稀疏性為後**，非窮舉式掃描。每個維度的分箱都在 v2 §2 明訂過界線，本研究不重新搜尋閾值。

### D1. LOH_Subtype（LOH 強度帶）

**實際來源**：由 C++ pipeline（`src/core/RegionProcessor.cpp` 的 `determine_loh_subtype()`）從兩個 per-region 信號組合而來：

**(a) `potential_loh` — HP 平衡度**（`SignificanceAnalyzer.cpp:316-323`）
```
hp_ratio = (n_HP1_family + 1) / (n_HP1_family + n_HP2_family + 2)   # Laplace smoothed
potential_loh = (hp_ratio < 0.1) OR (hp_ratio > 0.9)
```
即 HP1/HP2 family read 的比例極度偏頗（<10% 或 >90%）→ haplotype 失衡 → 疑似 LOH。

**(b) `verification_class` — Label-First × Cluster-First 雙軸交叉驗證**（`SignificanceAnalyzer.cpp:307-340`）

| Label-sig | Cluster-sig | verification_class |
|----------|-------------|-------------------|
| ✓ | ✓ | `Strong` |
| ✗ | ✓ | `Subclone` |
| ✓ | ✗ | `Weak` |
| ✗ | ✗ | `Noise` |

- Label-sig = 任一 label 維度（HP tag / Allele / Cluster）的 PERMANOVA 顯著
- Cluster-sig = 通過 gate 後 `global_alt.fisher p ≤ 0.05` 或 `global_hp_family.fisher p ≤ 0.05`

**合成真值表**：

| Potential_LOH | VerificationClass | → LOH_Subtype | 生物學意義 |
|---------------|-------------------|---------------|-----------|
| False | any | `None` | 無 LOH 訊號 |
| True | Noise | `LOH_Noise` | HP 失衡但 methylation 無顯著差異（訊號落 noise threshold 下） |
| True | Weak | `LOH_Weak` | HP 失衡 + label 有差但 cluster 結構弱 |
| True | Strong | `LOH_Strong` | HP 失衡 + 雙軸皆顯著（最純 LOH） |
| True | Subclone | `LOH_Subclone` | HP 失衡 + cluster 分群但 label 未對應（亞克隆 LOH） |

**重要修正**：v1 曾以 `LOH_Bed_Annotation` 為來源描述，實測該欄在當前 master TSV 中全為空；實際判定完全由 C++ per-region 的 HP_Ratio + VerificationClass 決定。本研究所用的 5 類都是 C++ pipeline 直出，並非外部 LOH bed overlay。

**為什麼保留 5 類而不合併**：合併 `Noise/Weak` 在 v2 §3 顯示會造成 LOH_Subtype 對 TP 判別力消失（粗 binning 吞掉訊號）。

**已知限制**：`Potential_LOH` 的 HP_Ratio 門檻為 self-phasing 輸出 — TO mode 的 HP tag 品質較差，可能導致 LOH_Subtype 在 TO 下有系統性偏移（見 memory `project_self_phasing_causal_chain_confirmed.md`）。

### D2. AF_class（Allele Frequency 粗類）

來源：`AF` 欄 → `af_class_coarse()` 函數。

| 類別 | 界線 |
|------|------|
| `Extreme` | AF < 0.1 OR AF > 0.9 |
| `Near-half` | 0.4 ≤ AF ≤ 0.6 |
| `Intermediate` | 其餘 |

**為什麼用三類而非 10-bin 細分**：10-bin 細切在 §11 實驗下顯示 n_cell 過低造成 Wilson CI 寬度 >0.2；三類提供足夠統計力同時保留「極端/中間/對等」三個生物意義區間。

### D3. cn_tier_F（SEQC2-grounded CN 分級）

來源：`CovM_used` 欄（diploid-normalized coverage multiple）→ `assign_cn_tier()` 用策略 F 的邊界 `[0.65, 0.99, 1.33, 1.82]`。

邊界來源：`kde_fix_validation/step1_hcc1395_seqc2/paired_pileup/per_cn_bin_metrics.tsv` 的 CN bin 中點。

| Tier | CovM_used 範圍 | 近似 CN |
|------|--------------|--------|
| T0 | <0.65 | sub-diploid / loss |
| T1 | 0.65–0.99 | CN=1–2 邊界 |
| T2 | 0.99–1.33 | CN=2 |
| T3 | 1.33–1.82 | CN=3 |
| T4 | ≥1.82 | CN≥4 |

### D4. HPFineNGroups（NG）

來源：`HPFineNGroups` 欄（C++ HPFineNGroups module 輸出，`--germline-hp-only` flag=off 模式下）。

類別：1, 2, 3, 4（4 為上限）。

**重要警告**（memory `project_hpfinengroups_subclone_marker.md`）：NG≥3 訊號在 flag=off 下由 somatic HP tag 人為分群驅動，flag=on 下會完全消失。本研究仍使用 flag=off 的值，因為 v2 既定 KDE baseline 都基於這個模式。

### D5. nr_band（NumReads 分箱）

來源：`NumReads` 欄 → `nr_band_label()` 三類。

| Band | NumReads 範圍 |
|------|-------------|
| `low` | <60 |
| `mid` | 60–119 |
| `high` | ≥120 |

**理由**：<60 覆蓋不足以穩定估計 ISM 指標；120 以上 KDE rescaling 已飽和（見 kde_fix_validation）。

---

## 2.2 Scheme 定義（S1-S5 白名單候選）

完整定義在 `ng_kde_rescaling/scripts/step6_tpfp_detailed.py` 的 `SCHEME_DEFS`。

| Scheme | 定義 | 生物學假設 |
|--------|------|-----------|
| **S1** `LOH_Strong_Extreme` | `LOH_Subtype == "LOH_Strong" AND AF_class == "Extreme"` | 純 haplotype loss，FP 極低 |
| **S2** `Subclonal_LOH_Inter` | `LOH_Subtype ∈ {LOH_Subclone, LOH_Weak} AND AF_class == "Intermediate"` | 亞克隆 LOH 中間 AF |
| **S3** `Diploid_Het` | `LOH_Subtype == "None" AND AF_class == "Near-half" AND cn_tier_F ∈ {T1, T2}` | CN=1-2 雜合體，典型 somatic het |
| **S4** `NonLOH_Extreme` | `LOH_Subtype == "None" AND AF_class == "Extreme"` | 無 LOH 極端 AF，germline 洩漏高風險桶 |
| **S5** `Combo` | `(S1 OR S2 OR S3) AND NOT S4` | 三個高 TP scheme 聯集，排除 S4 覆蓋 |

S6/S7 加 `HPFineNGroups>=3` 進一步限制，但在 HCC1395 TO 的 n 太小（`step6` 標 NO_USE），本研究不在主敘述中使用。

### 為什麼 S4 是關鍵診斷桶

S4（None + Extreme）在 HCC1395 TO 包含 21,652 TP + 8,780 FP = 30,432 variants（79.5% TO total），其 TP rate = 71.1% 正好等於 baseline。這就是「**無判別力區塊**」— 如果 scheme 比 baseline 沒差，代表這整個區塊的 biology（non-LOH 極端 AF）在 ISM 參數下沒有被切開。E5 的 sub-scheme 實驗就是試圖在 S4 內用單一 feature 再切，結果 fold≤1.17×（駁斥 H2）。

---

## 2.3 5D Cube 聯合空間

將 D1-D5 直接聯合，上限 `5 × 3 × 5 × 4 × 3 = 900` cells；實際活躍（n≥1）約 300-400 cells。

- **稀疏性控制**：Pareto/envelope 分析僅保留 `n >= 20` 的 cells（Wilson 95% CI 寬度仍可達 0.3，作為 weak evidence tier）
- **排序**：`tp_rate desc` 為主鍵，`n desc` 為 tie-break
- **累積封套（cumulative envelope）**：按排序疊加 cells 的 `n_tp`, `n_fp`，計算累積 purity、recall、fold-improvement

輸出位置：
- `data/tpfp_5d_cube_HCC1395_TO.tsv`（原始 cube）
- `data/tpfp_5d_pareto_HCC1395_TO.tsv`（排序後 ≥20 cells）
- `data/tpfp_5d_cumulative_envelope_HCC1395_TO.tsv`（累積封套）

---

## 2.4 統計指標

### 2.4.1 Precision / TP rate（純度）

\[ \text{TP rate} = \frac{n_\text{TP}}{n_\text{TP} + n_\text{FP}} \]

本研究**不計算 F1**，因為 master TSV 沒有 FN（caller-missed truth）資料。所有「判別力」指標用 `TP:FP ratio` 或 `purity` 表示。

### 2.4.2 TP:FP 比值與 fold-improvement

\[ \text{ratio} = \frac{n_\text{TP}}{n_\text{FP}}, \quad \text{fold} = \frac{\text{ratio}_\text{scheme}}{\text{ratio}_\text{baseline}} \]

若 `n_fp == 0`，ratio 記為 `inf`。在表中以 `>100×` 或 `∞` 顯示。

### 2.4.3 Recall

| 指標 | 計算 |
|------|-----|
| `tp_recall` | n_tp_in_scheme / n_tp_total（該 scheme 捕獲多少 TP） |
| `fp_recall` | n_fp_in_scheme / n_fp_total（該 scheme 捕獲多少 FP） |

TP recall 越高代表 scheme 覆蓋越廣；FP recall 越低代表 scheme 排除 FP 越有效。理想 scheme 是**高 TP recall、低 FP recall**。

### 2.4.4 Wilson 95% CI

`step6_tpfp_detailed.py` 使用 `scipy.stats` 二項 Wilson CI 估計 TP rate 的置信區間。`ci_width = ci_hi - ci_lo`，作為統計力指標。

### 2.4.5 Power tags（E3）

依 n 分級：

| 標籤 | n 範圍 | 用途 |
|------|-------|-----|
| `NO_USE` | <20 | 不可用，Wilson CI 過寬 |
| `WEAK` | 20-49 | 參考用，不納入主結論 |
| `OK` | 50-199 | 可用，但次要 |
| `STRONG` | ≥200 | 主結論引用 |

### 2.4.6 Mann-Whitney U + Cliff's delta（E4）

用於檢驗同一 scheme 內 TP vs FP 的 feature 分佈差異：

- `mannwhitney_p`: 兩側 U test p-value
- `cliffs_delta`: 取值 [-1, 1]，|δ|<0.1 無差，|δ|≥0.3 為中等效應

---

## 2.5 資料流程

```
master TSV (merged_7samples_paired_full_plus_hcc1395_to.tsv.gz)
  │
  ├─ step4_tpfp_discrimination.py      → S1-S7 on paired_full 7 samples
  ├─ step4b_tpfp_tomode.py             → S1-S7 on HCC1395 TO
  ├─ step6_tpfp_detailed.py            → E1-E5 metrics (fold, power, feature diffs, sub-schemes)
  └─ step8_multidim_panorama.py        → 5D cube + Pareto envelope
       │
       └─ figures/existing/fig_v2_1..11b
       └─ data/tpfp_*.tsv, tpfp_5d_*.tsv

obs01-obs08_*.py (本資料夾新增)
  ├─ obs01: AF 分佈 TP vs FP per sample
  ├─ obs02: CovM 分佈 TP vs FP
  ├─ obs03: 4-feature pairplot HCC1395 TO
  ├─ obs04: chr 空間分佈（S3 / Top-17 / Top-28 疊加）
  ├─ obs05: 單維度邊際 TP rate 貢獻
  ├─ obs06: NG 分佈 per sample 堆疊
  ├─ obs07: Pareto 軌跡細節
  └─ obs08: S3/S5/Top-17/Top-28 集合重疊
```

---

## 2.6 方法學限制

1. **S3 用 T1∪T2**：CN=1 理論上是 hemizygous；將其併入 S3 可能污染 Diploid Het 定義。v2 §3 驗證顯示 T1 的 TP rate 與 T2 接近（差異 <2%），合併合理。
2. **NG 訊號非 pure biology**：見 2.1 D4 警告。本研究使用 flag=off NG 值，所有結論需要 flag=on 下重驗證（out-of-scope）。
3. **5D cube 無平滑**：相鄰 cell 間可能出現 purity 抖動，empirical 排序 top-k 有過擬合風險（見 `06_limitations_future_work.md` §2）。
4. **Baseline TP:FP 依 caller 輸出**：不是真 F1 baseline，是 caller 判定的 TP/FP 比值。此研究測量的是「**條件過濾後 caller 輸出的 TP:FP 改進**」，不等於「過濾後真實 precision 改進」（因為 caller 本身可能漏報 TP）。
