<!--
建立時間: 2026-04-23 03:50
更新時間: 2026-04-23 03:50
狀態: in_progress
作者: InterSubMod Research
資料來源:
  - /big7_disk/liaoyoyo2001/InterSubMod/research/ng_kde_rescaling/data/merged_7samples_paired_full_plus_hcc1395_to.tsv.gz
  - scripts/analysis/20260423_B7_loh_noise_signal.py
相關假說: B7（Thread B Tier 2 觀察 §5.5：LOH_Noise + Extreme AF TP 88-96%）
驗證層級: L3（跨樣本重複觀察 + CI + Cohen's h）
Tier: T2（觀察確認；受限於「Extreme 由 HP_Ratio 結構決定」circular）
-->

# B7：LOH_Noise + Extreme AF 意外訊號的機制

## 1. 背景與假設

Thread B 週報觀察：`LOH_Noise` 子集在 Extreme AF 上的 TP rate 達 88–96%，與 `LOH_Strong` 幾乎持平。若此觀察為真，可能意味 `LOH.bed` 對「Noise-class LOH」的信心標記過於保守，該區域其實是 real LOH，只是訊號特徵被 `VerificationClass` 誤判。

本分析目標：

1. 在 HCC1395 TO（主軌 baseline TP=0.711）上確認 `LOH_Noise + Extreme` TP rate
2. 拆 CN tier × AF_class 驗證是否為混淆
3. 擴展到 7 樣本 paired_full + 可用 TO，檢驗跨樣本一致性
4. 結論 `LOH_Noise` 是否併入 `LOH_Strong` / 維持獨立 / 樣本特定

## 2. 資料與方法

- 主合併表：`merged_7samples_paired_full_plus_hcc1395_to.tsv.gz`（7 樣本 paired_full + HCC1395 TO/to_pileup + HCC1395_DORADO/H1437/H2009/HCC1937/HCC1954 TO）
- `LOH_Subtype` 定義（`RegionProcessor.cpp:188-204`）：
  - `Potential_LOH = HP_Ratio ≥ 0.85`（BAM HP tags）
  - `LOH_Noise = Potential_LOH ∧ VerificationClass=="Noise"`
  - `LOH_Strong = Potential_LOH ∧ VerificationClass=="Strong"`
- `AF_class`：Extreme / Intermediate / Near-half（VCF AF 切分）
- `CN_tier`：從 `Coverage_Category` lumped 為 `CN_le1 / CN_eq2 / CN_eq3 / CN_ge4`
- 分析腳本：`scripts/analysis/20260423_B7_loh_noise_signal.py`

## 3. 關鍵結果

### 3.1 HCC1395 TO：LOH_Noise ≈ LOH_Strong @ Extreme AF（意外）

表 1：HCC1395 TO 主軌 `Extreme AF` 下 pooled-CN TP rate（baseline TP=0.711）

| LOH_Subtype  | n      | k      | TP rate | Wilson 95% CI        |
|--------------|-------:|-------:|:-------:|:---------------------|
| LOH_Noise    |    431 |    401 | **0.930** | [0.902, 0.951]   |
| LOH_Strong   |    292 |    263 | **0.901** | [0.861, 0.930]   |
| LOH_Subclone |     24 |     23 |  0.958  | [0.798, 0.993]      |
| LOH_Weak     |  1,254 |  1,113 |  0.888  | [0.869, 0.904]      |
| None (非 LOH)| 30,432 | 21,652 |  0.711  | [0.706, 0.717]      |

**三件事同時成立**：
- Noise 與 Strong 的 CI 高度重疊（0.902–0.951 vs 0.861–0.930）
- 兩者相對 non-LOH baseline 都大幅提升（+0.19 ~ +0.22 TP rate）
- Cohen's h (Noise − Strong) = **+0.107**（small effect，不顯著）

### 3.2 LOH_Noise 幾乎全部落在 Extreme AF（關鍵結構限制）

表 2：LOH_Noise 的 AF_class 分布（7 樣本 × 可用 modes）

| sample × mode                    | Extreme | Intermediate | Near-half |
|----------------------------------|--------:|-------------:|----------:|
| COLO829 paired_full              |   1,483 |            0 |         0 |
| H1437 paired_full                |   1,968 |            0 |         0 |
| H1437 to_pileup                  |  15,251 |            0 |         0 |
| H2009 paired_full                |   7,138 |            0 |         0 |
| H2009 to_pileup                  |  41,267 |           20 |         0 |
| HCC1395 paired_full              |   4,637 |            1 |         0 |
| HCC1395 to_pileup (**主軌**)     |     431 |            0 |         0 |
| HCC1395_DORADO paired_full       |   4,027 |            0 |         0 |
| HCC1395_DORADO to_pileup         |  12,532 |            2 |         0 |
| HCC1937 paired_full              |   2,734 |            0 |         0 |
| HCC1937 to_pileup                |  10,392 |            0 |         0 |
| HCC1954 paired_full              |   1,022 |            0 |         0 |
| HCC1954 to_pileup                |   8,593 |            2 |         0 |

**解讀**：`Potential_LOH = HP_Ratio ≥ 0.85` 在 BAM 層面必然對應高 AF（homozygous 或 near-hom），因此 `LOH_Noise`（以及 `LOH_Strong/Weak`）基本上 **定義上就幾乎全部落在 Extreme AF**。這不是「Extreme AF 的特殊性讓 Noise 類 TP 變高」，而是「Extreme AF 是 LOH_* 子集的 default bucket」。

**推論**：B7 觀察 "LOH_Noise + Extreme TP 88-96%" 在 HCC1395 TO 成立（0.930），但與「LOH_Noise 整體 TP 0.930（Extreme 佔 100%）」同義，並非 Extreme 挑選出特殊 subset。

### 3.3 CN tier 拆分：HCC1395 TO 模式下 Noise 與 Strong 幾乎無差

表 3：HCC1395 TO `Extreme AF` × CN_tier × LOH_Subtype TP rate

| LOH_Subtype | CN_tier | n   | TP rate | Wilson CI         |
|-------------|---------|----:|:-------:|:------------------|
| LOH_Noise   | CN_le1  |  25 | 0.800 | [0.609, 0.911]      |
| LOH_Noise   | CN_eq2  | 229 | 0.930 | [0.890, 0.957]      |
| LOH_Noise   | CN_eq3  | 171 | 0.947 | [0.903, 0.972]      |
| LOH_Noise   | CN_ge4  |   6 | 1.000 | [0.610, 1.000]      |
| LOH_Strong  | CN_le1  |   8 | 0.750 | [0.409, 0.929]      |
| LOH_Strong  | CN_eq2  | 114 | 0.956 | [0.901, 0.981]      |
| LOH_Strong  | CN_eq3  | 159 | 0.893 | [0.835, 0.932]      |
| LOH_Strong  | CN_ge4  |  11 | 0.545 | [0.280, 0.787]      |

**解讀**：在每個 CN tier 內，Noise 與 Strong 的 CI 大幅重疊；HCC1395 TO 的 LOH_Noise Extreme 並非 CN-driven artifact。**Noise 在 HCC1395 TO 確實像 Strong 一樣是 high-TP 區域**。

### 3.4 跨樣本一致性：結論高度樣本特定

表 4：Noise vs Strong TP rate @ Extreme AF，per sample × mode

| sample             | mode         | Noise TP (n)        | Strong TP (n)      | Cohen's h (N−S) |
|--------------------|--------------|---------------------|--------------------|----------------:|
| COLO829            | paired_full  | 0.904 (n=1,483)     | 0.902 (n=244)      | **+0.007**      |
| H1437              | paired_full  | 0.999 (n=1,968)     | 0.999 (n=1,182)    | −0.006          |
| H2009              | paired_full  | 0.999 (n=7,138)     | 1.000 (n=3,516)    | −0.037          |
| HCC1395            | paired_full  | 0.970 (n=4,637)     | 0.997 (n=347)      | −0.238          |
| HCC1395_DORADO     | paired_full  | 0.992 (n=4,027)     | 1.000 (n=500)      | −0.174          |
| HCC1937            | paired_full  | 0.997 (n=2,734)     | 0.985 (n=546)      | +0.128          |
| HCC1954            | paired_full  | 0.992 (n=1,022)     | 1.000 (n=273)      | −0.175          |
| HCC1395            | to_pileup    | **0.930** (n=431)   | 0.901 (n=292)      | **+0.107**      |
| H1437              | to_pileup    | 0.747 (n=15,251)    | 0.877 (n=3,803)    | −0.336          |
| H2009              | to_pileup    | 0.909 (n=41,267)    | 0.955 (n=1,842)    | −0.188          |
| HCC1395_DORADO     | to_pileup    | 0.643 (n=12,532)    | 0.875 (n=3,364)    | −0.560          |
| HCC1937            | to_pileup    | **0.431** (n=10,392) | 0.822 (n=984)      | **−0.840**      |
| HCC1954            | to_pileup    | **0.251** (n=8,593)  | 0.394 (n=1,489)    | −0.306          |

**兩種行為模式**：

1. **Paired_full**：7/7 樣本 Noise ≈ Strong（|h|≤0.24，多為 negligible/small），且 TP rate 均 ≥0.90。**paired_full 下 Noise tier 確實被低估，可視為 real LOH 的 caller-confirmed 子集**。

2. **TO (pileup)**：樣本分裂嚴重：
   - 與 Strong 接近（|h|≤0.20）：HCC1395 (0.930), H2009 (0.909) — 2/6
   - 與 Strong 差距大（h<−0.30）：HCC1937 (0.431), HCC1954 (0.251), HCC1395_DORADO (0.643), H1437 (0.747) — 4/6
   - **TO 下 Noise 類 TP rate 跨樣本標準差巨大（0.251 → 0.930，range=0.679），中位數僅 0.691**。

## 4. 機制解讀

### 4.1 為何 HCC1395 TO 的 Noise 仍像 Strong？

- HCC1395 為**高純度 tumor-only 主軌**，在 ClairS-TO PASS 變體集合上，variant caller 已把大部分 FP 過濾；進到 LOH_Noise 的 region 本身的 `HP_Ratio≥0.85` 已是強 allele imbalance 訊號
- `VerificationClass="Noise"` 的 BAM-層訊號分類（entropy / permanova / ASM 組合）**偏向 read-level 的甲基化訊號弱**，但 **variant-level 不失準** — HCC1395 的 LOH.bed 本身可信（SEQC2 Jaccard=0.928）
- 因此 HCC1395 TO 下 Noise ≡ 「variant level real LOH + 甲基化不典型」，不是 「不信任的 LOH」

### 4.2 為何 HCC1937/HCC1954 TO 下 Noise 崩潰？

- HCC1937 TO FP rate 57%、HCC1954 TO FP rate 75% — 這兩個樣本本身 TO caller FP rate 異常高（相比 HCC1395 的 30%）
- Noise 區域在這些樣本中大量混入 non-LOH artifact（可能為 CNV_Loss / amplicon artifact 造成 HP_Ratio 偽高），VerificationClass 因無 consistent 甲基化訊號而落入 Noise
- Strong 類仍要求 verification 通過 → Strong 區仍可信但數量大降

### 4.3 與既有結論一致性

- 與 [07_LOH_CN_AF 研究總整理] 結論一致：`LOH` filter 方向已全面關閉（10/10 策略失敗），但 LOH Characterization 有效
- 與 Thread D NG=2 LOH-constrained phasing 觀察不衝突（LOH_Noise 與 NGroups 獨立）
- 與 `canonical_filter: NG=4 + AF<0.4 + NR≥80` 正交（此 filter 針對 FP rescue，非 LOH tier）

## 5. 對 "LOH.bed 是否過於保守" 的回答

**是，也不是**：

- **Paired mode**：LOH_Noise 其實絕大部分是 real LOH（TP rate ≥0.97），與 Strong 等價 → 「Noise」標籤在 paired 下確實低估。但 **此事對 filter 無用**，因為 paired Noise 已在高 TP 區（過濾反而損 TP）
- **TO mode**：4/6 樣本 Noise 崩潰（TP<0.75），Noise 確實包含真實不可信區 → **Noise 不能併入 Strong**，必須保留三 tier 區分
- **HCC1395 TO 單樣本是特例**（高純度 + 強 LOH.bed Jaccard），外推性有限

## 6. 決策：Noise tier 處理方式

| 路徑                                    | 建議   | 理由 |
|-----------------------------------------|--------|------|
| **併入 S1 擴展（全域）**                | ❌ NO-GO | TO 下 4/6 樣本 Noise TP<0.75，會引入大量 FP |
| 維持 Noise 獨立 tier                    | ✅ 採用 | paired 下無害、TO 下必須分離 |
| HCC1395-local「Noise 當 Strong」        | ⚠️ 不建議 | 單樣本 pilot、外推不穩、與既有 canonical filter 無顯著增益 |
| **Zone-Aware characterization annotation** | ✅ 採用 | 報告中標示 "paired_full 下 Noise≡Strong"，TO 下逐樣本驗證 |

**主結論**：B7 觀察 `LOH_Noise + Extreme AF TP 88-96%` **在 paired_full 全 7 樣本成立（TP≥0.90）**，但該「意外訊號」本質是 `LOH_*` 子集在定義上幾乎全部落在 Extreme AF（HP_Ratio≥0.85 ⇒ 高 AF），而非 Extreme AF 從 Noise 中挑出特殊 subset。TO mode 下 Noise 行為高度樣本特定，**不建議將 Noise 併入 S1**，維持現行 4-tier 分類。

## 7. 產出與限制

- **產出**：
  - `scripts/analysis/20260423_B7_loh_noise_signal.py`
  - `docs/experiments/in_progress/2026/04/figures/20260423_B7_loh_noise/fig_b7_hcc1395_to_heatmap.png`
  - `docs/experiments/in_progress/2026/04/figures/20260423_B7_loh_noise/fig_b7_cross_sample_extreme.png`
  - `research/ng_kde_rescaling/data/b7_loh_noise_stratified_fullsamples.tsv`
  - `research/ng_kde_rescaling/data/b7_loh_noise_hcc1395_TO.tsv`
  - `research/ng_kde_rescaling/data/b7_loh_noise_per_sample.tsv`
  - `research/ng_kde_rescaling/data/b7_loh_noise_summary.json`

- **限制**：
  - HCC1395 TO LOH_Noise n=431 只是該樣本全量 40,115 的 1.1%，結論為區域特性而非 filter-level
  - TO pileup 模式的 Noise 崩潰可能與 `VerificationClass` 在 TO 下的 QS 失效相關（AUC=0.497）— 屬既有結論
  - 未區分 `LOH_Bed_Overlap`（LOH.bed 命中）子集內的 Noise，可進一步深挖但非當前 blocker
