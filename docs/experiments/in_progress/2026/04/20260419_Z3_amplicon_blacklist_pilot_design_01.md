---
title: Z3 × HCC1954 Amplicon/High-CN Blacklist Pilot — 設計文件
date: 2026-04-19
status: design
owner: InterSubMod Research
scope: TO mode · HCC1954-specific FP · Method A（Amplicon / High-CN Blacklist）
related:
  - 20260418_Z3_Internal_Feature_Exploration_01.md
  - ../../../../research/z3_internal_feature_exploration/00_INDEX.md
---

# Z3 × HCC1954 Amplicon/High-CN Blacklist Pilot — 設計

## 一、背景與問題定義

### 1.1 Z3 Zone 定位

Zone-Aware Confidence Framework 將 TO mode call 依 `loh_flag × af_extreme × HPFineNGroups` 分為 Z1–Z5。**Z3** 定義為：

```
Z3 = (TO LOH hit) ∧ (caller_af < 0.1 ∨ caller_af > 0.9) ∧ (HPFineNGroups ≤ 1)
```

跨樣本 TP rate σ=0.28（H2009 0.868 → HCC1954 0.050），是 TO mode 中最不穩定的 Zone（詳見 Z3 Internal Feature Exploration 研究，20260418 報告）。

### 1.2 Z3 內部無跨樣本特徵可用

先前 Step 1/2.5 研究結論 **NEGATIVE**：

- 12 候選特徵 × 7 樣本 AUC 掃描，無特徵在 ≥3 樣本達 |AUC|≥0.60
- AF∈[0.4,0.6] × CN × NGroups germline band 分層僅 1/7 樣本符合 pattern

→ Z3 內部無「全局可用」的二階分層特徵。

### 1.3 HCC1954 的特殊 FP pattern（Step 3 觀察）

HCC1954 Z3 FP **93% 集中於 chr5 (43%) / chr8 (29%) / chr17 (13%)**，配合：

- FP NumReads=55 vs TP=37（Mann-Whitney p=4.7e-9）
- FP Coverage_Multiple=0.73 vs TP=0.49（CN proxy ≈ 1.5×）

![Z3 chr 分佈對比（HCC1954 vs HCC1395）](../../../../research/z3_internal_feature_exploration/figures/z3_chr_distribution_contrast.png)

這對應 HCC1954 已知的 HER2+/pseudo-tetraploid breast cancer 架構（ERBB2 focal amp + arm-level 8p loss + 17p TP53 LOH）。

---

## 二、假說與方法選擇

### 2.1 本 pilot 假說

**H-ABL**：HCC1954 Z3 FP 大宗來自高 CN / amplicon / arm-level LOH 區的 read mis-phasing / AF artifact；在 Z3 內部用「基因體位置黑名單」或「sample-specific High-CN 閾值」過濾，可在**不破壞其他 6 樣本**的前提下拉高 HCC1954 F1。

### 2.2 為什麼選 Method A（Blacklist）而非其他

| 方法 | 放棄原因 |
|------|---------|
| B. 甲基化層甲基化-only 重分類 | Z3 NG≤1 → 甲基化特徵結構性退化（Step 1 驗證） |
| C. 全樣本 AF 0.4–0.6 germline rule | H-Z3d 僅 1/7 樣本符合（Step 2.5 NEGATIVE） |
| D. per-sample ML | 樣本數太少、對論文故事無幫助 |
| **A. Amplicon/High-CN Blacklist** | **最直接、可解釋、可用文獻座標驗證合理性** |

### 2.3 合理性檢核（Sanity Check）

#### (a) Z3-conditional 限制（**避免破壞 non-Z3 TP**）

- 若對全 call 套用 blacklist，會把 chr5/8/17 上高信心 non-Z3 TP 全部誤殺
- 本 pilot 所有 strategy 以 `mask & in_z3 == True` 只對 Z3 內部啟用

#### (b) 文獻 focal coord vs 實際 arm-level FP 分佈（Step 4 校核）

將文獻常報的 HER2/MYC focal coord 套入 HCC1954 Z3 FP 檢驗：

| 文獻座標 | FP 捕獲 | 實際 FP 分佈 |
|---------|--------|-----------|
| chr17q12 HER2 (35–42Mb) | **0/583** | chr17 FP 53% 在 chr17p (0–20Mb, TP53 LOH 區) |
| chr8q24 MYC (120–140Mb) | 37/583 | chr8 FP 68% 在 chr8p (0–40Mb, 8p loss 區) |
| chr5p gain (1–50Mb) | 132/875 | chr5 FP 分佈全長 |

→ **實際 FP 是 arm-level 架構問題，非單純 focal amplicon**。因此 S1 blacklist 改為 arm-level（TP53 LOH、8p loss、MYC arm、HER2 arm、整段 chr5）而非文獻狹窄焦點。

![HCC1954 Z3 FP Pos 分佈 vs 文獻 amplicon](../../../../research/z3_internal_feature_exploration/figures/step4_hcc1954_fp_pos_hist.png)

#### (c) Circularity Guard（S3 only）

S3 若用 HCC1954 自己的 Coverage_Multiple 分佈定 cutoff → 等於用 outcome 定 filter，循環定義。

**修正策略**：cutoff 用**非 Z3 region 的 baseline 95th percentile** 計算（與 FP/TP 無關），再套回 Z3 過濾。

#### (d) Ceiling Calibration（S4）

加入「拒絕全 Z3」的 ceiling strategy，計算**最大可能 gain**。任何 refined strategy 若低於 ceiling 太多 → 無競爭力；若接近 ceiling → refined 的 collateral damage 才有意義。

---

## 三、Strategy 定義

### S1 — Literature arm-level ∩ Z3

基於 HCC1954 已知 cytogenetic 架構：

```python
S1_LITERATURE_ARMS = [
    ("chr5",           1_000_000, 180_000_000, "chr5 (HCC1954 複合重排)"),
    ("chr8",                   0,  45_000_000, "chr8p loss"),
    ("chr8",         120_000_000, 145_000_000, "chr8q MYC amplicon"),
    ("chr17",                  0,  25_000_000, "chr17p TP53 LOH"),
    ("chr17",         35_000_000,  42_000_000, "chr17q12 HER2 focal"),
]
```

### S2 — Whole chr5/8/17 ∩ Z3

對 chr5、chr8、chr17 整條染色體上所有 Z3 call 全部拒絕。代表「極端粗粒度」，用於評估「是否需要精細 arm 切分」。

### S3 — Sample-specific High-CN ∩ Z3

每個樣本以其自身**非 Z3 region** 的 `Coverage_Multiple` 95th percentile 為 cutoff：

```
cutoff[sample] = quantile(CovM[sample, ~in_z3], 0.95)
blacklist = (CovM > cutoff[sample]) & in_z3
```

- 優點：sample-specific，不需文獻座標
- 風險：HCC1954 Z3 本身 CovM 偏高 → cutoff 偏高 → 捕獲有限（預期 TP 極少數，FP 也少）

### S4 — Reject all Z3（Ceiling）

直接把 HCC1954 的全 Z3 region 標為 FN，計算 F1 上限。

---

## 四、F1 計算規則

```
TP_kept = TP ∧ ¬mask
FP_kept = FP ∧ ¬mask
TP_lost = TP ∧ mask
FN_new  = FN_original + TP_lost

Precision = TP_kept / (TP_kept + FP_kept)
Recall    = TP_kept / (TP_kept + FN_new)
F1_new    = 2·P·R / (P+R)
ΔF1       = F1_new − F1_before
```

其中 `FN_original = 0`（ClairS-TO 對 regions-of-interest 報出的 call 無獨立 missed truth 列），故 `Recall` 退化為 `TP_kept / (TP_kept + TP_lost) = 1 − TP_lost / TP_total`。

---

## 五、預期結果矩陣

| Strategy | HCC1954 期望 ΔF1 | 其他 6 樣本期望 | 合理性 |
|---------|----------------|--------------|------|
| S1 arm-level ∩ Z3 | 正（接近 ceiling） | 輕微負（chr5/8/17 少量 Z3 正確 call 誤殺） | 高（文獻支撐） |
| S2 whole-chr ∩ Z3 | 正（略高於 S1） | 輕微負（範圍更大） | 中（過度粗） |
| S3 CovM 95%ile ∩ Z3 | 小正 | 近零 | 高（無循環） |
| S4 reject all Z3 | 最大正（ceiling） | 大負 | 僅作 ceiling，非可採用 |

**成功判準**：

- **POSITIVE**：S1 或 S3 在 HCC1954 ΔF1 ≥ +0.005 且其他 6 樣本 mean ΔF1 ≥ −0.003
- **CONDITIONAL**：HCC1954 有增益但其他樣本 mean ΔF1 ∈ [−0.005, −0.003]
- **NEGATIVE**：HCC1954 增益 < +0.002 或 collateral damage mean ΔF1 < −0.005

---

## 六、產物規劃

| 類型 | 路徑 | 用途 |
|------|------|------|
| 腳本 | `research/z3_internal_feature_exploration/scripts/step4_pos_distribution_check.py` | 文獻座標 vs 實際 FP 位置校核 |
| 腳本 | `research/z3_internal_feature_exploration/scripts/step5_amplicon_blacklist_pilot.py` | S1–S4 F1 計算 |
| 數據 | `data/step4_hcc1954_fp_amplicon_hit.tsv` | 文獻座標 hit rate |
| 數據 | `data/step5_blacklist_pilot_results.tsv` | per-sample × strategy ΔF1 |
| 數據 | `data/step5_blacklist_pilot_summary.tsv` | 匯總表 |
| 圖表 | `figures/step4_hcc1954_fp_pos_hist.png` | FP 位置 vs 文獻 amplicon |
| 圖表 | `figures/step5_blacklist_delta_f1.png` | ΔF1 per-sample × strategy |
| 報告 | `20260419_Z3_amplicon_blacklist_pilot_result_01.md` | 結果 + 合理性驗證 + 最終判定 |

---

## 七、風險與限制

1. **HCC1954 Z3 FP=50218，TP=17068**：FP 絕對數字巨大。即使 blacklist 僅捕獲 4% FP，仍可能使 precision 明顯上升
2. **其他 6 樣本的 collateral damage**：chr5/8/17 不僅 HCC1954 有 CNV，其他樣本也可能有局部 amp（如 H2009 EGFR、COLO829 CDKN2A 等），大範圍 blacklist 會誤殺
3. **Zone-Aware Framework 定位不變**：即使本 pilot POSITIVE，仍屬 HCC1954-specific conditional rule，非跨樣本通用 filter
4. **與 F pilot 的關係**：F pilot（`NG=4+AF<0.4+NR≥80`）已覆蓋 HCC1954 大部分 CNV artifact。本 pilot 的 Z3-conditional blacklist 是 **補充性**而非取代
5. **合理性的最後檢驗**：S1/S2 捕獲的 TP（~71–77 個）必須 spot-check 是否**真 TP**（SEQC2 truth set 內），若為 truth set annotation 錯誤 → 實際 ΔF1 更高
