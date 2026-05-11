---
title: 敘述比較分析 — 從 Baseline 到 5D Pareto 的層層逼近
status: in_progress
last_updated: 2026-04-22
---

# 04. 敘述比較分析（主文）

> **核心結論**：在 HCC1395 TO mode（baseline TP=71.1%）下，從 biology-informed scheme (S3=95.5%, fold 8.69×) 到聯合 5D 空間 Top-17 envelope (96.1%, fold **10.0×**)，存在明顯的 Pareto 改進；S3/S5 既非唯一高純度點，也不是 Pareto-optimal。但此改進尚未跨樣本驗證（H4 未證）。

---

## 4.1 故事線概覽

```
Layer 0: Baseline
  HCC1395 TO baseline TP:FP = 2.47:1 (TP rate 71.1%)
  COLO829 paired_full baseline TP:FP = 15.5:1 (TP rate 93.9%)
        │
        ▼
Layer 1: 單維度切分
  LOH_Subtype 單獨？ AF_class 單獨？ CN_tier 單獨？
  → obs05 顯示 Near-half AF 單獨提 TP rate 從 71% → ~93% (但 recall 很低)
        │
        ▼
Layer 2: Biology-Informed Scheme (v2 §3)
  S3 Diploid Het (None + Near-half + T1/T2) = 95.5% TP
  S5 Combo (S1|S2|S3) & !S4 = 91.8% TP
  → fold 8.69× (S3) 或 4.53× (S5)
        │
        ▼
Layer 3: 單 feature threshold 在 scheme 內再切 (v2 §14, E5)
  最佳 S4d = S4 + 單特徵 threshold → fold 最多 1.17× (over S4 baseline)
  → 🛑 scheme 內 refinement 飽和 (H2 駁斥)
        │
        ▼
Layer 4: 5D 聯合空間 (v2 §17, E6)
  Top-17 cells = 96.1% purity + 3.7% recall (fold 10.0×)
  Top-28 cells = 90.4% purity + 7.4% recall (fold 4.73×)
  → ✅ Pareto-dominate S3 (purity 更高 + recall 2.9× 更多)
        │
        ▼
Layer 5: 限制與未驗證
  H4 跨樣本泛化未驗證；COLO829 top cells 座標結構不同 (NG=1 vs NG=2-3)
```

---

## 4.2 Layer 0: Baseline 的意義

Caller 自己判定的 TP:FP 比值是本研究的對照基準：

| 測試場 | Baseline TP:FP | Baseline TP% | 備註 |
|-------|--------------:|-----:|------|
| HCC1395 TO | **2.47:1** | 71.1% | 唯一有足夠 FP (n=11,606) 做判別力量測；TO mode 特有的 germline 洩漏 |
| COLO829 paired_full | 15.5:1 | 93.9% | 次要測試場；FP=2,273 可做驗證 |
| HCC1395 paired_full | 47.5:1 | 97.9% | FP 稀少，頂板效應 |
| 其他 5 樣本 paired_full | 64-8433:1 | >98.5% | 不可量測改進 |

**為什麼 HCC1395 TO 是主測試場**：
1. FP 基數最大（11,606）— 統計力足以區分 fold 1.5× 和 8×
2. Baseline TP% 只有 71.1% — 留下 29% 改進空間
3. TO mode 代表實際使用 ClairS-TO 的 FP 生成機制（無 normal 參考的 somatic caller）

---

## 4.3 Layer 1: 單維度切分（obs05）

將 HCC1395 TO 按單一維度切分，對比 baseline 的 TP rate 變化：

| 維度 | 最高 TP rate 類別 | 該類別 n | 該類別 TP rate | 相對 baseline (71.1%) 變化 |
|------|------------------|---------|---------------|--------------------------|
| AF_class | Near-half | 937 | **93.9%** | **+22.8%** |
| LOH_Subtype | LOH_Strong | 1,087 | 89.6% | +18.5% |
| cn_tier_F | T2 | 11,783 | 78.2% | +7.1% |
| HPFineNGroups | NG=2 | 13,226 | 75.3% | +4.2% |
| nr_band | mid | 22,318 | 73.9% | +2.8% |

**解讀**：
- **AF_class 是最強單維度判別器**（Near-half 提 23 個百分點）
- **LOH_Subtype 次之**（Strong/Weak 都顯著高於 None）
- CN_tier 和 NG/NR 的單維度影響較弱（<10%）
- 沒有任何單維度能把 TP rate 推過 94%；需要**組合**

![obs05 單維度邊際 TP rate 貢獻](figures/new/obs05_dim_contribution.png)

![obs01 per sample AF 分佈（TP vs FP 疊加）](figures/new/obs01_af_distribution_per_sample.png)

![obs02 per sample CovM 分佈](figures/new/obs02_cn_distribution_per_sample.png)

**這直接暗示**：若只用一個維度做 filter，最多只能到 ~94% purity，無法超越 biology-informed 組合（>95%）。

---

## 4.4 Layer 2: Biology-Informed Scheme（v2 §3 + §11）

| Scheme | n | TP rate | TP:FP ratio | fold vs baseline | TP recall |
|--------|--:|-----:|----------:|----------------:|---------:|
| **S3** Diploid_Het | 380 | **95.5%** | 21.35 | **8.69×** | 1.3% |
| **S1** LOH_Strong_Extreme | 292 | 90.1% | 9.07 | 3.69× | 0.9% |
| **S2** Subclonal_LOH_Inter | 214 | 87.4% | 6.93 | 2.82× | 0.7% |
| **S5** Combo (S1|S2|S3) & !S4 | 886 | **91.8%** | 11.14 | **4.53×** | 2.9% |
| S4 NonLOH_Extreme | 30,432 | 71.1% | 2.47 | 1.00× | 76% |

**S3 的三個條件**：`LOH_Subtype == "None" AND AF_class == "Near-half" AND cn_tier_F ∈ {T1, T2}`

**S4 的診斷意義**：S4 占 TO mode 79.5% 的 variants，TP rate 恰好 = baseline；這個「無判別力桶」是 E5 sub-scheme 試圖再切的目標。

![v2_3 filter scheme bar（S1-S8 TP rate / recall / FP reduction）](figures/existing/fig_v2_3_filter_scheme_bar.png)

![v2_5 scheme operating points（purity vs recall）](figures/existing/fig_v2_5_operating_points.png)

![v2_7 baseline vs scheme TP:FP fold-improvement](figures/existing/fig_v2_7_baseline_tp_fp_fold.png)

![v2_4 paired_full 7 樣本 × S1-S6](figures/existing/fig_v2_4_per_sample_schemes.png)

### 4.4.1 S3 vs S5 的 trade-off

S3 是最乾淨的 operating point（95.5% 純度），但只捕 1.3% TP。S5 犧牲 3.7 個百分點換 2.3× 的 recall。如果目標是「高純度白名單」，選 S3；如果目標是「中等純度但較廣覆蓋」，選 S5。

---

## 4.5 Layer 3: Sub-scheme Refinement 飽和（v2 §14, E5）

**假設 H2**：在 S4 內加單一 feature threshold 可切出 sub-high-TP。

**實驗**：在 S4 內分別加 NumReads≥80, NumReads<80, NG=2, NG=3, 等 9 個候選 threshold。結果：

| Best sub-scheme | 條件 | n | TP rate | fold vs S4 baseline |
|----------------|------|--:|------:|-------------------:|
| S4d | S4 + NG=2 + high NR | 3,242 | 76.7% | **1.17×** |
| S4a-S4c, S4e-S4i | 各種組合 | - | <76% | <1.17× |

**結論**：**H2 被駁斥**。S4 內單特徵 refinement 最大 fold 僅 1.17×，遠低於 S3 的 8.69×。這意味著 S4 這個 biology 條件（non-LOH + 極端 AF）本身就是 FP 富集區，無法用 feature threshold 挽救。

**但這不等於整個參數空間飽和**。E5 的失敗引導 E6：若放棄「在 scheme 內加 threshold」框架，改成**聯合 5D cube**，能否突破？

![v2_9 S4 內 TP vs FP feature violin（AlleleDelta δ=-0.49）](figures/existing/fig_v2_9_feature_violin_in_S4.png)

![v2_10 sub-scheme operating points（fold 飽和證據）](figures/existing/fig_v2_10_subscheme_operating_points.png)

---

## 4.6 Layer 4: 5D Pareto Envelope（v2 §17, E6）**主要突破**

將 `LOH_Subtype × AF_class × cn_tier_F × HPFineNGroups × nr_band` 直接聯合為 900 cells，保留 n≥20 的活躍 cell，按 TP rate 降序累積：

### 4.6.1 關鍵 operating points

| Operating Point | n | cum TP rate | cum TP recall | TP:FP ratio | fold vs baseline |
|-----------------|--:|----------:|------------:|------------:|----------------:|
| Baseline (HCC1395 TO) | 40,115 | 71.1% | 100% | 2.47 | 1.00× |
| S3 Diploid Het | 380 | 95.5% | 1.3% | 21.35 | 8.69× |
| **Top-17 cells** | **1,099** | **96.1%** | **3.7%** | **24.56** | **10.0×** ⭐ |
| S5 Combo | 886 | 91.8% | 2.9% | 11.14 | 4.53× |
| **Top-28 cells** | **2,285** | **90.4%** | **7.4%** | **11.68** | **4.73×** |
| Top-32 (probe) | 8,421 | 83.0% | 24.5% | 4.89 | 1.98× |

### 4.6.2 Top-17 Pareto-dominates S3

| 指標 | S3 | Top-17 | 勝負 |
|------|---:|-------:|----|
| Purity | 95.5% | **96.1%** | Top-17 勝 |
| Recall | 1.3% | **3.7%** | Top-17 勝 (2.8×) |
| TP:FP | 21.35 | **24.56** | Top-17 勝 |
| fold | 8.69× | **10.0×** | Top-17 勝 |

**兩個指標同時更好 = Pareto dominance**（下圖 Panel B，S3/S5 標記落在 Top-17 外側即代表 Pareto-dominated）。

![obs07 Pareto 軌跡：purity vs k + recall-purity envelope](figures/new/obs07_pareto_trajectory.png)

![v2_11 5D cube panorama HCC1395 TO](figures/existing/fig_v2_11_panorama_HCC1395_TO.png)

### 4.6.3 Top-17 內的 cell 組成（見 obs04 + obs08）

前 17 個 cell 的生物特徵（來自 `data/tpfp_5d_pareto_HCC1395_TO.tsv`）：

```
rank  LOH_Subtype    AF_class       cn_tier  NG  NR_band  n    TP%
 1    LOH_Strong     Extreme        T1       3   mid      22   100.0%
 2    LOH_Noise      Extreme        T3       2   high     56    98.2%
 3    LOH_Weak       Extreme        T2       3   mid      44    97.7%
 4    LOH_Weak       Intermediate   T1       2   mid      43    97.7%
 5    None           Near-half      T1       2   mid     171    97.7%  ← S3 subset
 6    LOH_Noise      Extreme        T2       2   high    146    97.3%
 7    LOH_Noise      Extreme        T2       3   high     35    97.1%
 8    None           Near-half      T2       3   high     31    96.8%  ← S3 subset
 9    LOH_Strong     Extreme        T2       3   high     87    96.6%
10    LOH_Weak       Intermediate   T1       3   mid      27    96.3%
...
```

**兩個令人驚訝的 pattern**：

1. **LOH_Noise + Extreme + high_NR** 出現 4 次（rank 2, 6, 7 + 後續）— 這是 S1/S2 都未捕獲的組合。`LOH_Noise` 傳統被視為 noise tier，但配合 high NR + 特定 NG 組合後 purity 97%+。
2. **LOH_Weak + Intermediate + T1** 連續出現（rank 4, 10）— S2 捕獲 LOH_Weak+Intermediate 但未細分 T1/T2；細分後 T1 subset 達 97.7% TP。

### 4.6.4 集合重疊分析（obs08）

來自 `tables/obs08_scheme_overlap_matrix.tsv`：

| | S3 | S5 | Top-17 | Top-28 |
|-|--:|--:|------:|------:|
| S3 size | **380** | - | 337 (89%) | 366 (96%) |
| S5 size | - | **886** | 540 (61%) | 730 (82%) |
| Top-17 size | - | - | **1,099** | 1,099 |
| Top-28 size | - | - | - | **2,285** |

**關鍵解讀**：
- **Top-17 有 1,099 regions，但只 337 (31%) 屬 S3** — 多出的 762 regions 是 S3 無法捕獲的新高純度區
- **Top-28 有 2,285 regions，只 730 (32%) 屬 S5** — 多出的 1,555 regions 是 S5 無法捕獲的
- 5D cube 不是把 S3/S5 「切更細」，而是**發現 S3/S5 根本沒碰到的高純度 cell**

![obs08 scheme 集合重疊 heatmap + size stack](figures/new/obs08_overlap_venn_schemes.png)

![obs04 Top-17 / Top-28 / S3 在 chromosomes 上的空間分佈](figures/new/obs04_chr_spatial_scheme.png)

![obs03 HCC1395 TO 4-feature pairplot（AF × CovM × NR × AlleleDelta）](figures/new/obs03_feature_pairplot.png)

---

## 4.7 跨樣本對照：COLO829 paired_full

在 COLO829 paired_full（baseline 93.9%，FP=2,273）重跑相同 5D cube：

| Operating Point | Purity | Recall | Fold vs baseline 15.5:1 |
|-----------------|-----:|------:|-----------------------:|
| COLO829 baseline | 93.9% | 100% | 1.00× |
| COLO829 S5 | 94.7% | 6.1% | 1.15× |
| COLO829 Top-k cells | >98% (n 小) | <5% | 1.3-1.5× |

**已知差異**（已在 v2 §17.3、§17.7 記錄）：
- COLO829 top cells 的主軸是 **NG=1**（而 HCC1395 TO 是 **NG=2-3**）
- COLO829 高純度 cell 多在 **Intermediate AF**（而 HCC1395 TO 是 **Near-half**）
- 這兩個樣本的「最優 cell 座標」結構上不同 → **H4 (top-k 跨樣本泛化) 未證**

### 4.7.1 為什麼座標不同

推測：
- COLO829 baseline 已 94%，多數 FP 來自 germline leak 在 **CN-neutral LOH + Intermediate AF** 區
- HCC1395 TO baseline 71%，FP 主要來自 **non-LOH Extreme AF**（S4 area）
- 所以兩個樣本的「FP 熱區」本來就不同，**不該期待 top cells 完全吻合**

但這也代表：**不能在 HCC1395 TO 上訓練 top-k list 後直接用於 COLO829**。跨樣本部署需要每個 sample 自己的 top-k（或跨樣本共識 cells）。

![v2_11b COLO829 paired_full 5D panorama（座標主軸與 HCC1395 TO 不同）](figures/existing/fig_v2_11b_panorama_COLO829.png)

![v2_8 paired_full 7 樣本 × scheme TP rate heatmap](figures/existing/fig_v2_8_per_sample_scheme_heatmap.png)

![obs06 per sample NG 分佈 TP vs FP（hatched）](figures/new/obs06_ng_by_sample_stacked.png)

---

## 4.8 圖表速查

| 圖 | 內容 | 支援哪層結論 |
|----|------|--------------|
| `figures/new/obs01_af_distribution_per_sample.png` | AF 分佈 TP vs FP | Layer 1 — 驗證 Near-half AF 偏 TP |
| `figures/new/obs02_cn_distribution_per_sample.png` | CovM 分佈 | Layer 1 — CN=2 (T2) 是最大桶 |
| `figures/new/obs03_feature_pairplot.png` | 4-feature pairplot | Layer 1 — 看二維投影是否可分 |
| `figures/new/obs04_chr_spatial_scheme.png` | chr 空間 S3/Top-17/Top-28 分佈 | Layer 4 — 驗證非集中在單染色體 |
| `figures/new/obs05_dim_contribution.png` | 單維度邊際 TP rate | Layer 1 — 單維最強為 AF |
| `figures/new/obs06_ng_by_sample_stacked.png` | NG 堆疊 | Layer 4 — NG 分佈跨樣本差異 |
| `figures/new/obs07_pareto_trajectory.png` | Pareto 軌跡 Top-17/28/32 | Layer 4 — Top-17 P-dominates S3 |
| `figures/new/obs08_overlap_venn_schemes.png` | S3/S5/Top-17/Top-28 重疊 | Layer 4 — 5D 新發現 ~500 regions |

既有 v2 圖（`figures/existing/`）：
- fig_v2_3 filter_scheme_bar — Layer 2 scheme 總覽
- fig_v2_4 per_sample_schemes — Layer 2 跨樣本 S1-S7
- fig_v2_5 operating_points — Layer 2 purity-recall 散點
- fig_v2_7 baseline_tp_fp_fold — Layer 2 fold-improvement
- fig_v2_9 feature_violin_in_S4 — Layer 3 TP/FP 特徵差
- fig_v2_10 subscheme_operating_points — Layer 3 sub-scheme 飽和證據
- fig_v2_11 panorama — Layer 4 5D cube 全景

---

## 4.9 主結論（五點；2026-04-22 rerun 後更新 H4 判決）

1. **H1 強支持**：LOH × AF × CN 能切分 TP/FP（S3 達 95.5%，fold 8.69×）
2. **H2 駁斥**：在預定 scheme 內加單特徵 threshold 無法明顯再提 purity（最佳 S4d fold 1.17×）
3. **H3 強支持**：5D 聯合空間存在 Pareto-optimal 點 — Top-17 cells 同時在 purity (96.1% vs 95.5%) **與** recall (3.7% vs 1.3%) 兩指標上勝過 S3
4. **Top-17 發現 ~760 regions S3/S5 捕不到的新高純度 cell** — 5D 不是 refinement，是 **enlargement**
5. **H4 升級為中-強支持（原未驗證）**：

### 4.9.1 H4 升級證據表（2026-04-22 ISM rerun 後新增）

| 證據 | rerun 前（HCC1395 TO vs COLO829 pf）| rerun 後（6 TO samples 對比）| H4 判決 |
|------|-----------------------------------|--------------------------|---------|
| LOH gradient 跨樣本 | HCC1395-only | ✅ **6/6 TO 樣本方向一致** | ↑ 升級 |
| NG=3 boost 跨樣本 | HCC1395-only | ✅ **5/6 TO 樣本可見** | ↑ 升級 |
| 2D heatmap 色梯度 | 僅 HCC1395_TO | ✅ **5/6 TO 樣本**（HCC1954 反向例外）| ↑ 升級 |
| 3D cells 跨樣本 `n_samples_high ≥ 5` | 未計算 | ✅ **16 cells**（Hard Gate 閾值 10，**超標 60%**）| ↑ 完全觸發 |
| Top-17 cells 在 16-cell 共識集 | 未計算 | ⚠ **0/17 重疊** | Top-17 仍為 HCC1395-specific |
| 跨樣本共識主軸 | 未知 | `Extreme AF + NG=2 + T1-T3 CN`（非 Top-17 的 Near-half + NG=3）| 新發現 |

### 4.9.2 16-cell 共識集建議（取代 Top-17 作為 LOSO 目標）

新研究重點從 HCC1395 Top-17 轉移至**跨樣本 16-cell 共識集**：
- Tier A（3 cells）：`LOH_Strong/Weak/None + Extreme + T1 + NG=2/4 + mid_NR`，達 7/13 樣本
- Tier B（6 cells）：`Extreme + T0-T4 + NG=2-4`，達 6/13
- Tier C（7 cells）：全部 NG=2，達 5/13

詳見 `07_figure_layers.md` §7.5.1 Tier A/B/C 表與 §7.7 Hard Gate。

**LOSO 設計建議**：
- 保留原 HCC1395 Top-17 作為 baseline（內部最優）
- 新增 16-cell consensus set 作為 primary filter（外部泛化性）
- 比較兩者在 held-out 樣本上的 F1 gap（預期 Top-17 崩塌 >50% / 16-cell 保留 >70%）

**跨文件引用**：本結論的原始數據與完整推論見 `07_figure_layers.md` §7.5-7.7；TSV 數據 `data/cell_consistency_matrix.tsv`。

---

## 4.10 為什麼這篇分析比 v2 報告更有價值

v2 報告（20260422_02.md）已涵蓋技術細節，但本資料夾做了三件 v2 沒做的事：

| v2 報告 | 本研究增補 |
|---------|----------|
| fig_v2_11 panorama 1 張圖 | **obs07 Pareto 軌跡細節** — 放大 Top-17/28/32 的三個 operating point 比較 |
| 未量化 S3 vs Top-17 重疊 | **obs08 + overlap matrix TSV** — 明確算出 Top-17 多出 762 regions (69%) 不在 S3 |
| 未檢視 chr 空間分佈 | **obs04** — 確認 Top-17 cells 分散於所有 chromosomes，非單點集中 |
| 單維度討論零散 | **obs05** — 量化 5 維度各自的邊際 TP rate 貢獻 |
| 跨樣本對照僅有 fig_v2_11b | **obs06 NG 堆疊** — 視覺呈現 7 樣本 × 2 mode 的 NG 分佈差異 |
| 特徵差異僅有 S4 內 violin | **obs03 pairplot** — 全 TO mode 4-feature 聯合投影 |

本研究的定位是**將 v2 的結論實體化為可重現、可重跑的分析資料夾**，並加強可視化以支持後續 reviewer 或 PI 閱讀。

---

## 4.11 下一步

1. **立即**（本資料夾範疇）：
   - 完成 `05_biology_interpretation.md`（top-17 cells 的生物機制詮釋）
   - 完成 `06_limitations_future_work.md`（H4 未證 + 論文定位）
2. **短期**（1-2 週內）：
   - LOSO 檢驗：在 HCC1395 TO 找到的 top-k cells 應用到 COLO829，測是否保留 fold >1.2×
   - 若 LOSO 失敗 → 研究跨樣本共識 top cells（intersection of sample-specific top lists）
3. **中期**（依 PI 決策）：
   - 若 reviewer 要求 F1 → 取得 SEQC2 FN 資料重算
   - 若 ML 方向升級 → logistic regression 在 5D 座標上，替代 empirical top-k
