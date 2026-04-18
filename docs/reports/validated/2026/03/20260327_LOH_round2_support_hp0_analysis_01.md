<!--
建立時間: 2026-03-27
勘誤: 2026-03-30 — TO 端全面修正（HP integer tag bug）：Core1 tier enrichment、Core2 HP0 by LOH status、§5.2 HP0 比較、§5.3 HP0 by tier（C0: 9.7%→85.3%）、§5.4 Stage 4 affinity（neutral: 95%→75%）、§6.1-6.3 結論更新、§7 意外發現更新。見各節 [修正] 標注
目標: LOH Round 2 正式報告 — effective_hp support 分層 + HP0 來源分析
處理範圍: 748,391 region rows, 7 樣本, paired + TO mode
關聯檔案:
  - docs/reports/validated/2026/03/20260327_LOH_round1_cross_sample_audit_01.md
  - docs/architecture/20260327_InterSubMod研究願景定錨_01.md
  - docs/plans/2026/03/20260327_LOH_round2_execution_spec_01.md
  - output/synthesis/observation_workspaces/20260327_loh_round2_support_hp0_analysis/
  - scripts/analysis/build_loh_round2_support_hp0_analysis.py
-->

# LOH Round 2 正式報告：effective_hp Support 分層 + HP0 來源分析

> 生成時間：2026-03-27
> 性質：observation-first，結論含多個意外發現（unexpected findings）
> 分析腳本：`scripts/analysis/build_loh_round2_support_hp0_analysis.py`
> 輸出 workspace：`output/synthesis/observation_workspaces/20260327_loh_round2_support_hp0_analysis/`

---

## 一、這輪要回答的問題

Round 1 確立了「LOH-like 分佈與 TP/FP 差異診斷底圖」，但遺留了：

1. **Core 1**：排除 `effective_hp_reads < 10` 的雜訊後，Tier A（≥30 reads）的 LOH-like FP/TP enrichment 是否更穩定、更具泛化力？
2. **Core 2**：TO 中大量 HP0 的成因是什麼？HP0 在 LOH-like region 是否特別高？Stage 4 affinity score 是否有補充 evidence 的潛力？

---

## 二、Tier 定義與資料範圍

### 2.1 Tier 定義（本輪固定）

| Tier | `effective_hp_reads` | 說明 |
|------|---------------------|------|
| C0 | = 0 | 完全無 HP reads，`HP_Ratio=0.5` 是假象 |
| C  | 1–9 | 極弱 support，幾乎全是雜訊（本輪確認） |
| B  | 10–29 | 中等 support，需小心解讀 |
| A  | ≥ 30 | 強 support，可信 LOH evidence 候選 |

### 2.2 資料範圍

- 資料來源：Round 1 `all_region_rows.tsv.gz`（748,391 rows）
- 有效 TP/FP rows：748,391（TP=616,580, FP=131,811）
- 7 個樣本（paired + TO 各一組）

---

## 三、Core 1 結果：effective_hp Support 分層

### 3.1 全域 Tier 分析 — Paired mode

| Tier | TP Total | FP Total | TP LOH-like% | FP LOH-like% | FP/TP Enrichment | Fisher p | 解釋 |
|------|---------|---------|------------|------------|----------------|---------|------|
| C0 | 29,873 | 135 | 0.0% | 0.0% | — | 1.0 | 完全無 HP signal |
| C | 2,296 | 50 | 67.5% | 70.0% | 1.037× | 0.76 | 統計無顯著差異 |
| **B** | 25,073 | 1,299 | **38.6%** | **34.8%** | **0.901×** | **0.006** | ⚠️ 反轉！FP LOH-like 反而比 TP 低 |
| **A** | 268,028 | 1,945 | **31.4%** | **36.7%** | **1.169×** | **7.2e-07** | ✅ 顯著 FP enrichment |

> **Round 1 整體 paired enrichment = 1.194×**，對比顯示：**Tier A 幾乎是整體 enrichment 的全部來源。**

### 3.2 全域 Tier 分析 — TO mode

~~舊版表格（HP bug，已廢棄）~~：
~~C0: 27,016 / 12,783~~  ~~C: 45,914 / 18,092 (0.990×)~~  ~~B: 72,903 / 29,558 (1.002×)~~  ~~A: 145,477 / 67,949 (0.867×)~~

> **[修正 2026-03-30]** TO 端 eff_hp 全部重算，tier 分組數量大幅改變：

| Tier | TP Total | FP Total | TP LOH-like% | FP LOH-like% | FP/TP Enrichment | Fisher p | 解釋 |
|------|---------|---------|------------|------------|----------------|---------|------|
| C0 | 732 | 338 | 0.0% | 0.0% | — | 1.0 | 完全無 HP signal |
| C（<10）| 4,084 | 2,027 | 76.1% | 77.3% | 1.015× | 0.32 | 低 phasing，微弱無差異 |
| B（10-29）| 27,422 | 14,506 | 55.8% | 52.4% | **0.938×** | — | TP 微弱富集 |
| **A（30-49）** | **44,177** | **25,255** | **61.4%** | **43.3%** | **0.706×** | — | **TP 富集最強** |
| **A+（≥50）** | **214,895** | **86,256** | **39.1%** | **30.0%** | **0.766×** | — | **TP 富集** |
| A≥30（合計）| 259,072 | 111,511 | 42.9% | 33.0% | **0.769×** | p≈0 | TP 富集，極顯著 |

> **修正後整體 TO enrichment = 0.805×**（TP LOH 44.5% vs FP LOH 35.8%）。TO Tier A 和 A+ 均呈 TP 富集（<1），方向一致。

### 3.3 Paired Tier A per-sample

| Sample | TP Total | FP Total | TP LOH-like% | FP LOH-like% | FP/TP Enrichment | Fisher p | 說明 |
|--------|---------|---------|------------|------------|----------------|---------|------|
| HCC1954 | 16,695 | 27 | 9.4% | 33.3% | 3.531× | 5.8e-04 | 強訊號但 FP 數量少 |
| H2009 | 127,961 | 85 | 28.7% | 76.5% | 2.664× | 1.1e-19 | 最可靠的高 enrichment |
| H1437 | 45,522 | 4 | 29.6% | 75.0% | 2.534× | 0.08 | **FP 數量不足（n=4），不可信** |
| HCC1937 | 12,040 | 189 | 55.5% | 84.7% | 1.526× | 2.8e-17 | 中等 enrichment，顯著 |
| HCC1395_DORADO | 24,174 | 198 | 48.3% | 59.1% | 1.222× | 2.7e-03 | 輕度 enrichment |
| COLO829 | 14,630 | 900 | 9.9% | 11.8% | 1.188× | 0.076 | **不顯著** |
| HCC1395_5kHz | 27,006 | 542 | 46.5% | 46.9% | 1.007× | 0.90 | **完全無效果** |

### 3.4 TO Tier A（30-49）per-sample

~~舊版表格（HP bug，全部無效）~~

> **[修正 2026-03-30]** 修正後 TO Tier A per-sample 結果完全不同（COLO829 從 232→10,881 TP regions，H1437 從 11,144→8,635）：

| Sample | TP Total | FP Total | TP LOH-like% | FP LOH-like% | FP/TP Enrichment | Fisher p |
|--------|---------|---------|------------|------------|----------------|---------|
| HCC1395 5kHz | 7,277 | 2,867 | 84.4% | 82.9% | **0.982×** | 0.062 |
| HCC1395_DORADO | 6,784 | 2,680 | 86.7% | 84.9% | **0.979×** | 0.022 |
| **COLO829** | **10,881** | **6,066** | **15.6%** | **14.8%** | **0.951×** | 0.190（不顯著）|
| H1437 | 8,635 | 2,586 | 79.8% | 75.4% | **0.945×** | 3.0e-6 |
| H2009 | 6,431 | 545 | 76.6% | 74.5% | **0.973×** | 0.270（不顯著）|
| HCC1937 | 645 | 629 | 91.3% | 85.2% | **0.933×** | 8.8e-4 |
| HCC1954 | 3,524 | 9,882 | 28.0% | 25.3% | **0.904×** | 0.002 |

> **修正後結論**：TO Tier A 的 enrichment **全部 < 1（TP 富集）**，方向一致，與全域一致，矛盾不復存在。
> 舊版「HCC1395_DORADO 和 H1437 enrichment > 1.17×」完全是 HP bug 造成的假象（大量 HP:i:11/21 reads 未計入 HP1/HP2Family，導致 eff_hp 被低估，Tier A 分組成員有誤）。

---

## 四、Core 1 意外發現：Tier B 反轉是 Sample Composition Artifact

### 分析

Paired Tier B 整體 enrichment = 0.901×（顯著反轉）。但 per-sample 分析揭示：

| Sample | TP (Tier B) | FP (Tier B) | Enrichment |
|--------|------------|------------|-----------|
| **COLO829** | 18,088 | **1,206 (93%)** | 1.093× |
| HCC1395 | 1,944 | 56 | 0.841× |
| HCC1395_DORADO | 1,667 | 27 | 0.867× |
| H1437 | 782 | 1 | — |
| H2009 | 1,667 | 1 | — |
| HCC1937 | 220 | 6 | — |
| HCC1954 | 705 | 2 | — |

COLO829 貢獻了 93% 的 Tier B FP，且 COLO829 本身 Tier B enrichment = **1.093×（正向）**。

整體 Tier B 的反轉原因：
- 非 COLO829 樣本（HCC1395, HCC1395_DORADO 等）TP LOH-like fraction 很高（70-77%）
- 這些樣本幾乎沒有 FP，大量高 LOH-like TP 拉高了整體 TP LOH% (38.6%)
- 但 FP LOH% (34.8%) 由 COLO829（93% FP）決定，維持在 32.9%

**結論：Tier B 整體反轉不是真實的生物訊號，是不同樣本之間的 TP/FP 比例差異導致的聚合 artifact。**

---

## 五、Core 2 結果：HP0 來源分析

### 5.1 LOH-like vs non-LOH-like 的 HP0 分佈

#### Paired mode（「真實 LOH」的參照）

| 分組 | HP0 mean | HP0>0 fraction | HP0>10% fraction |
|------|---------|--------------|----------------|
| TP / non-LOH-like | 8.97% | 46.3% | 22.1% |
| **TP / LOH-like** | **4.11%** | **26.3%** | **10.9%** |
| FP / non-LOH-like | 8.87% | 48.8% | 25.0% |
| **FP / LOH-like** | **6.60%** | **32.6%** | **17.0%** |

**Paired LOH-like → 顯著較低的 HP0**（TP: 4.1% vs 9.0%，p<1e-200）。
這與預期一致：若一個 region 真的是 one-sided phasing（LOH），大多數 reads 都在 HP1-family，留下 HP0 的就少。

#### TO mode（「雜訊型 LOH」的觀察）

~~舊版表格（LOH-like 分組依賴舊版 Potential_LOH，n 有誤）~~

> **[修正 2026-03-30]** LOH-like 分組（Potential_LOH）在修正後改變，n 減少（TP LOH-like: 175,630→129,589），但 HP0 數值本身（NHP0/NumReads）不受 HP bug 影響：

| 分組 | n | HP0 mean（修正後）| 方向 |
|------|---|----------------|------|
| TP / non-LOH-like | 161,721 | **4.56%** | — |
| **TP / LOH-like** | **129,589** | **9.57%** | ↑ LOH-like HP0 更高 |
| FP / non-LOH-like | 82,429 | **4.66%** | — |
| **FP / LOH-like** | **45,953** | **10.48%** | ↑ LOH-like HP0 更高 |

**🚨 TO LOH-like → HP0 仍然更高（修正後確認）**：TP（9.6% vs 4.6%），FP（10.5% vs 4.7%）

**HP0 數值可信（不受 HP bug 影響），趨勢方向不變。** 但數值從舊版（7.6%/5.5%）升高至修正後（9.6%/4.6%），差距更明顯。

### 5.2 解釋：Paired vs TO 的 HP0 行為為何相反？

| 模式 | LOH-like region 的 HP0 | 解釋 |
|------|----------------------|------|
| Paired | **低 HP0**（4.1%） | Real LOH：one-sided phasing 把大多 reads 指向一個 family，HP0 自然少 |
| TO | **高 HP0**（**9.6%**） [修正 2026-03-30] | Partial phasing artifact：TO phasing 在高 tumor purity 下形成的 LOH-like 訊號，伴隨 phasing 不完全 → 部分 reads 逃到 HP0 |

**這提供了一個潛在的 TO LOH-like 品質判斷指標**：
- TO LOH-like + 低 HP0 → 可能是較「clean」的 LOH 訊號
- TO LOH-like + 高 HP0 → 可能是 phasing artifact 產生的假 LOH

### 5.3 HP0 的 Tier 依賴性（重要品質確認）

| Mode | Tier | TP HP0 mean | FP HP0 mean | 說明 |
|------|------|------------|------------|------|
| paired | C0 | 40.6% | 28.7% | 完全靠 HP0 構成，無意義 |
| paired | **C** | **40.4%** | **54.1%** | **HP0 極高！Tier C LOH 完全不可信** |
| paired | B | 10.1% | 7.8% | HP0 仍高，LOH 訊號也不穩定 |
| paired | **A** | **3.3%** | **5.6%** | 可接受的 HP0 水平 |
| to | C0 | **85.3%** | **85.4%** | 幾乎全是 HP0，無 HP phasing（[修正 2026-03-30]） |
| to | **C** | **77.7%** | **77.8%** | **極高 HP0，Tier C LOH 完全不可信（[修正 2026-03-30]）** |
| to | B | **22.1%** | **20.9%** | 中高 HP0（[修正 2026-03-30]） |
| to | **A** | **9.2%** | **7.2%** | 可接受的 HP0 水平（[修正 2026-03-30]） |
| to | **A+** | **2.7%** | **2.2%** | 最低，phase 支持品質最好（[修正 2026-03-30]） |

> **[修正 2026-03-30]**：TO 行數值因 HP integer tag bug 修正而大幅更新。修正前 C0/C/B/A 各 Tier 的 HP0 值被嚴重低估（C0 舊值 9.7% → 修正後 85.3%，因為 HP:i:11/21 reads 過去未計入 eff_hp，使低 eff_hp tier 中混入大量本應在高 tier 的 regions）。修正後 HP0 單調遞減更合理（tier 越高 → eff_hp 越高 → HP0 越低）。

**Paired Tier C 的 HP0 高達 40-54%，TO Tier C 更高達 78%** — 兩者均顯示：在 1-9 有效 HP reads 的 region，大量 reads 是 unphased。任何 Tier C 的 LOH-like 訊號幾乎全是統計噪聲，應直接排除。修正後 TO Tier C0 更達 85% HP0，直接確認 Tier C0 在 TO mode 完全不可信。

### 5.4 Stage 4 Affinity（HP0 reads 傾向哪個 HP family？）

#### 全量（含 HP0=0 的 region）

幾乎所有 region 的 affinity score = 0（94-98% 是「中性」），因為大多數 region 根本沒有 HP0 reads。整體分析不具意義，需過濾 NHP0>0。

#### NHP0 > 0 的 region（Paired: 39.8% = 133,109 rows；TO: 39.2% = 164,461 rows）

| Mode | Truth | LOH-like | n | Affinity abs mean | Frac neutral | NHP0 median |
|------|-------|---------|---|-------------------|------------|------------|
| paired | TP | non-LOH | 106,511 | 0.070 | 90.9% | 7.0 |
| paired | TP | LOH | 25,120 | **0.054** | **90.6%** | 4.0 |
| paired | FP | non-LOH | 1,087 | 0.125 | 78.7% | 4.0 |
| paired | FP | LOH | 391 | **0.068** | **88.5%** | 5.0 |
| to | TP | non-LOH | **59,527** | **0.095** | **1.6%** | 4.0 |[修正 2026-03-30] |
| to | TP | LOH | **58,125** | **0.038** | **74.8%** | 7.0 |[修正 2026-03-30] |
| to | FP | non-LOH | **27,280** | **0.102** | **1.5%** | 4.0 |[修正 2026-03-30] |
| to | FP | LOH | **19,529** | **0.043** | **73.0%** | 8.0 |[修正 2026-03-30] |

> **[修正 2026-03-30]**：TO 行因 HP integer tag bug 修正而更新。Potential_LOH 分組成員改變（原本錯誤的 LOH 分組排除，正確 LOH 分組加入），affinity 計算也因 HP1/HP2Family 參考組更正而不同。

**關鍵觀察（TO，修正後）**：
- TO LOH-like 的 HP0 affinity 仍然較低（0.038-0.043），但只有 **73-75% neutral**（不是舊版誤報的 95%）
- TO non-LOH-like 的 HP0 affinity 明顯更高（0.095-0.102），且 **neutral 只有 1.5-1.6%**（修正後 non-LOH affinity 方向性更清晰）

**解釋**：在 TO LOH-like region 中，HP0 reads 方向性較弱（只有 25-27% 有方向），代表部分 phasing ambiguity。TO non-LOH-like region 的 HP0 reads 多數有方向性（98% 非中性），代表這些 region 的 reads 能被有效地分配到 HP1/HP2 family。這與「TO LOH-like 是 partial phasing 的產物」假說一致，但 Stage 4 的分辨力比原本預期略強。

如果 HP0 reads 在 LOH-like region 有強 affinity，可以補充 LOH evidence；但 TO LOH-like 的 affinity（25-27% 非中性）仍有限，**Stage 4 affinity 在 TO LOH-like region 的補充價值仍然有限**，不如 paired non-LOH FP 的 affinity（21% 非中性）。

### 5.5 H1 假說：HP0 ↔ Coverage 相關性

| Mode | Truth | Spearman ρ | p-value | 解釋 |
|------|-------|-----------|---------|------|
| paired | FP | **+0.247** | 7.7e-49 | **中等正相關** |
| paired | TP | +0.053 | 1.3e-202 | 弱正相關 |
| to | FP | -0.034 | 2.2e-34 | 極弱負相關 |
| to | TP | -0.034 | 6.2e-76 | 極弱負相關 |

**Paired FP 的特殊發現**：Coverage 越高，HP0 越多（ρ=0.247）。
這可能代表 paired FP region 傾向出現在「reads 多但 phasing 困難」的複雜區域（例如重複序列、CNV 邊界），而不是在高 coverage 就能解決 phasing 的正常區域。

### 5.6 H3 假說：同位點 TO HP0 > Paired HP0

| Sample | Median HP0 Delta (TO-paired) | Wilcoxon p | 結論 |
|--------|------------------------------|-----------|------|
| COLO829 | +0.000 | p≈0 | TO HP0 系統性高於 paired |
| H1437 | +0.000 | 1.7e-272 | TO HP0 系統性高於 paired |
| H2009 | +0.000 | 4.3e-44 | TO HP0 系統性高於 paired |
| HCC1954 | +0.000 | 4.5e-14 | TO HP0 系統性高於 paired |
| HCC1937 | +0.000 | 0.38 | **不顯著** |
| HCC1395 | +0.000 | 1.00 | **HCC1395 TO 的 HP0 與 paired 相當** |
| HCC1395_DORADO | +0.000 | 1.00 | **HCC1395_DORADO 亦相當** |

中位數 delta 均為 0（因為分佈是零膨脹的），但 Wilcoxon 在多個樣本顯示 TO HP0 tail 系統性高於 paired。

**HCC1395（兩版本）是例外**：TO 的 HP0 與 paired 沒有系統性差異。這可能代表 HCC1395 的 longphase-to phasing 品質相對較好，或者其 LOH landscape 在 TO 中沒有造成額外的 phasing 困難。

---

## 六、整合推論與結論

### 6.1 Core 1 主要結論

**結論 1**：Paired 端，LOH-like FP enrichment 完全集中在 Tier A（≥30 reads）。Tier B 的「反轉」是 COLO829 主導的聚合 artifact，不是真實的生物訊號。

**結論 2**：Paired Tier A 的 enrichment（1.169×）在統計上顯著，但 **per-sample 異質性仍然很大**（HCC1954 3.5× vs HCC1395 1.0×）。作為通用 filter 的能力有限，需 sample-aware 設計。

**結論 3**：TO 端，~~Tier A enrichment 的方向是 **sample-specific 的**。HCC1395_DORADO 和 H1437 在 TO Tier A 中 LOH-like FP 確實多於 TP（1.18–1.32×）~~。**[廢棄：此結論基於 HP bug 錯誤數據]**。**修正後（2026-03-30）**：ALL 7 個樣本的 TO Tier A（30-49 reads）enrichment 均 < 1（TP 富集），方向一致，COLO829 Tier A（0.951×, p=0.189）不顯著但方向相同；全域 TO Tier A 的 enrichment = **0.706×**（TP 富集），Tier A+（≥50）= **0.766×**（TP 富集）。TO 與 paired 方向**相反**（paired Tier A+ 是 FP 富集）。

**結論 4**：Tier C（1–9 reads）的 LOH-like 幾乎全是雜訊，Paired Tier C HP0 高達 40-54%，應在後續所有分析中直接排除。

### 6.2 Core 2 主要結論（HP0）

**結論 5（Paired vs TO 的根本差異）**：

| 觀察 | Paired | TO |
|------|--------|-----|
| LOH-like 的 HP0 | **低**（4.1%） | **高**（**9.6%**）[修正 2026-03-30] |
| 方向 | LOH-like = clean one-sided | LOH-like = partial phasing artifact |
| Stage 4 affinity（有HP0時）| 輕度非中性（LOH: 90% neutral） | 部分中性（LOH: **75% neutral**）[修正 2026-03-30] |

Paired LOH-like 是「real one-sided phasing」的反映。TO LOH-like 有相當比例是 phasing 不完全產生的，伴隨更高的 HP0 和無方向的 affinity。

**結論 6（Tier C HP0 量化確認）**：Paired Tier C 有 40-54% HP0；**TO Tier C 有 78% HP0，TO Tier C0 有 85% HP0**（修正後 2026-03-30，比舊版 15-19% 更高）。任何 Tier C 的 LOH 分析都應被標記為「不可信」，TO Tier C0/C 的結論比 Paired 更強。

**結論 7（Stage 4 affinity 的局限）**：Stage 4 affinity 在 TO LOH-like 區域方向性有限（**75% neutral**，修正後不是 95%），補充 LOH evidence 的能力仍然有限。修正後發現 TO non-LOH-like 的 HP0 affinity 極強（只有 1.5% neutral），這是新發現，但主要是 non-LOH 特性而非 LOH 相關，對 LOH filter 設計影響有限。只有 Paired FP region 的 HP0 affinity 有潛在利用價值（21% 非中性）。

**結論 8（樣本特異性）**：HCC1395（兩版本）在 H3 假說中不顯著（TO HP0 ≈ paired HP0），代表 HCC1395 的 phasing 在 TO 和 paired 之間品質相當。其他樣本（COLO829、H1437、H2009、HCC1954）則 TO HP0 系統性高於 paired。

### 6.3 LOH 的可用性更新（對照研究願景）

| 情境 | 可用性評估 | 條件 |
|------|-----------|------|
| Paired Tier A LOH-like → FP risk feature | **有條件可用** | 僅 HCC1954、H2009、HCC1937 有足夠 enrichment 且顯著；COLO829 和 HCC1395 無效 |
| TO Tier A LOH-like → TP/FP 分離 | **TP 富集（全樣本一致）**[修正 2026-03-30] | 修正後所有樣本方向一致：TO Tier A(30-49) 0.706×，Tier A+(≥50) 0.766×，與 paired 方向相反 |
| Tier C 任何 LOH 訊號 | **不可用** | Paired HP0 高達 54%，訊號全是雜訊 |
| TO LOH-like + 低 HP0 → 較可信 | **潛在可用，待 Round 3 驗證** | 低 HP0 TO LOH-like 可能是較 clean 的訊號 |
| Stage 4 affinity 補充 | **TO LOH-like 無用，Paired FP 有潛力** | Paired FP 的 HP0 affinity 相對較強（0.125） |

---

## 七、意外發現摘要

| 發現 | 預期 | 實際 | 重要性 |
|------|------|------|--------|
| TO Tier A enrichment 方向 | LOH-like 在 TO FP 較多 | **修正後（2026-03-30）全樣本一致 TP 富集（Tier A: 0.706×，Tier A+: 0.766×）；舊值（0.867×）是 HP bug 假象** | 🔴 高 — TO LOH 是 TP 指標而非 FP 指標（與 paired 相反） |
| Paired LOH-like HP0 vs 非 LOH-like | 相當 | **LOH-like HP0 低 50% 以上** | 🟡 中 — 驗證了 paired LOH 的真實性 |
| TO LOH-like HP0 vs 非 LOH-like | 與 paired 同方向 | **TO LOH-like HP0 反而更高** | 🔴 高 — 揭示 TO LOH-like 的部分 artifact 成分 |
| Stage 4 affinity 在 TO LOH-like | 有方向性，可補充 | **95% 中性，幾乎無用** | 🟡 中 — Stage 4 在 TO LOH 補充有限 |
| Tier C HP0 水平 | 較高但可接受 | **Paired Tier C: 40-54%!** | 🔴 高 — 直接排除 Tier C LOH 的任何解讀 |
| Tier B 反轉 | 表示中間 support 的 LOH 不穩定 | **Sample composition artifact（COLO829 主導）** | 🟢 低 — 非真實發現，是聚合方法論的問題 |

---

## 八、Round 2 的限制

1. **Tier 閾值敏感性**：本輪固定 Tier A=30, B=10，未做 A=20 或 A=50 的敏感性分析
2. **HP0 分析是觀察性的**：無法直接從資料因果說明 TO HP0 是 phasing failure 造成的
3. **Stage 4 affinity 需更多樣本**：目前結論基於 7 個樣本，需要更多 paired/TO 比較
4. **HCC1954 的 TO FP 偏高**：HCC1954 有 35,074 TO FP（其他樣本 TO FP 最多 17,528），主導了 TO 全域統計，這個不平衡可能影響解讀

---

## 九、下一步建議

### 9.1 直接可行的後續分析

1. **TO LOH-like 品質分層**：TO LOH-like + `hp0_ratio < 0.05`（低 HP0）vs `hp0_ratio > 0.10`（高 HP0）的 FP/TP enrichment 是否有差異？這是 Round 3 的主要問題
2. **Paired Tier A LOH-like FP filter 試驗**：對 H2009、HCC1937 這兩個有強 enrichment 且數量足夠的樣本，測試 `Tier A LOH-like AND hp0_ratio < 0.05` 的 filter 效果
3. **Tier 閾值敏感性分析**：用 Tier A=20 和 A=40 重跑 Core 1，確認結論的穩健性

### 9.2 與研究願景的對接

- **目標 1**（per-CpG 甲基多標籤）：LOH Tier 定義應作為 per-CpG 標籤的一個維度（Tier A LOH-like 才是有效的 HP 相關標籤）
- **目標 5**（F1 提升）：Paired Tier A LOH-like 對 H2009、HCC1937 可嘗試作為 FP risk feature，但需注意 COLO829 和 HCC1395 的無效性
- **目標 4**（TO normal 補強）：TO LOH-like + HP0 的組合特徵值得進一步探索，特別是「低 HP0 的 TO LOH-like」是否更接近 paired 的 clean LOH

---

## 十、圖表索引

**Fig01 — 全域 Tier FP/TP Enrichment（paired vs TO，Tier A/B/C/C0 對比）**

![Fig01 Tier Enrichment Global](../../../../../research/loh_investigation/figures/loh_round2/fig01_tier_enrichment_global.png)

**Fig02 — Paired Tier × Sample Enrichment Heatmap（各樣本 × Tier 的 paired FP/TP enrichment）**

![Fig02 Paired Tier Enrichment Heatmap](../../../../../research/loh_investigation/figures/loh_round2/fig02_paired_tier_enrichment_heatmap.png)

**Fig03 — Tier LOH-like Fraction Barplot（各 Tier 的 LOH-like fraction，TP vs FP 並排）**

![Fig03 Tier LOH-like Fraction Barplot](../../../../../research/loh_investigation/figures/loh_round2/fig03_tier_loh_frac_barplot.png)

**Fig04 — HP0 Distribution Violin（HP0 ratio 分佈，LOH-like vs non-LOH-like，paired vs TO）**

![Fig04 HP0 Distribution Violin](../../../../../research/loh_investigation/figures/loh_round2/fig04_hp0_distribution_violin.png)

**Fig05 — Stage 4 Affinity Score Distribution（HP0/HP3 reads 傾向 HP1 vs HP2 family 的分佈）**

![Fig05 Affinity Score Distribution](../../../../../research/loh_investigation/figures/loh_round2/fig05_affinity_score_distribution.png)

**Fig06 — Same-locus HP0 Delta（同位點 TO-paired HP0 差值分佈，by concordance group 與 sample）**

![Fig06 HP0 Same-locus Delta](../../../../../research/loh_investigation/figures/loh_round2/fig06_hp0_same_locus_delta.png)

**Fig07 — HP0 Ratio vs Coverage Scatter（H1 假說：低 coverage 是否導致高 HP0？）**

![Fig07 HP0 vs Coverage](../../../../../research/loh_investigation/figures/loh_round2/fig07_hp0_vs_coverage.png)

**Fig08 — HP0 Mean by Tier（各 Tier 的 HP0 平均值，TP vs FP，paired vs TO）**

![Fig08 HP0 Mean by Tier](../../../../../research/loh_investigation/figures/loh_round2/fig08_hp0_mean_by_tier.png)

---

## 十一、數據輸出索引

| 檔案 | 說明 |
|------|------|
| core1_tier_enrichment_global.tsv | 全域 mode × tier enrichment |
| core1_tier_enrichment_by_sample_mode.tsv | 各樣本 × mode × tier enrichment |
| core2_hp0_by_loh_status.tsv | HP0 統計 by LOH status |
| core2_hp0_by_sample_loh_status.tsv | HP0 統計 by sample × LOH status |
| core2_hp0_loh_mannwhitney.tsv | HP0 Mann-Whitney 結果 |
| core2_affinity_by_loh_status.tsv | Stage 4 全量 affinity 統計 |
| core2_affinity_nhp0gt0_only.tsv | Stage 4 affinity（NHP0>0 過濾後） |
| core2_affinity_loh_mannwhitney.tsv | Affinity Mann-Whitney 結果 |
| core2_hp0_coverage_correlation.tsv | HP0 vs coverage Spearman 相關 |
| core2_hp0_same_locus_delta.tsv | 同位點 HP0 delta 統計 |
| core2_hp0_by_concordance.tsv | HP0 by concordance group |
| core2_hp0_by_tier.tsv | HP0 分佈 by Tier |
