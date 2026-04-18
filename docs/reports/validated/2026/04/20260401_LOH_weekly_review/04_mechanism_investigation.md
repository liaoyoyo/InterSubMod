<!--
建立時間: 2026-04-01 12:00
目標: 週回顧第四章——TO vs Paired LOH enrichment 方向翻轉的機制調查與 TP rescue 潛力評估
處理範圍: 同位點 concordance 分析、HP balance 分佈、AF 中位數差異、VerificationClass 分層、低 AF TP rescue
關聯檔案:
  - docs/reports/validated/2026/03/20260330_post_hp_fix_to_loh_investigation_01.md
  - scripts/analysis/build_post_hp_fix_to_loh_investigation.py
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_post_hp_fix_to_loh_investigation/
-->

# 04. TO vs Paired LOH Enrichment 方向翻轉：機制調查

---

## 1. 核心問題

第 03 章確認了一個關鍵事實：HP 修正後 TO LOH enrichment = 0.805x（TP-enriched），而 Paired LOH enrichment = 1.194x（FP-enriched）。

**為什麼相同的基因組位點，使用不同的 calling mode（paired vs TO）會導致 LOH enrichment 方向完全翻轉？**

本章透過同位點 concordance 分析、HP balance 分佈、AF 分層等手段，逐步推導出機制解釋。

---

## 2. 同位點 LOH Concordance 分析（最關鍵證據）

### 2.1 方法

取 Paired 和 TO 都有 call 的**完全相同基因組位點**，比較兩種模式對同一位點的 LOH 判定結果。如果方向差異只是統計波動，同位點的 LOH 判定應大致一致；如果是系統性差異，應觀察到一致的偏向。

### 2.2 數據

總計 **288,609 個 valid 同位點 loci**（TP subset），涵蓋 7 個樣本。

#### Table 1. 同位點 LOH Concordance by Sample

| 樣本 | 同位點數 | TO_only_LOH | Paired_only_LOH | 比值 (TO/Paired) | Concordance |
|------|---:|---:|---:|:---:|:---:|
| H1437 | 44,843 | 5,802 | 112 | **51.8x** | 86.8% |
| COLO829 | 32,809 | 5,205 | 207 | **25.1x** | 83.5% |
| HCC1954 | 16,342 | 2,421 | 96 | **25.2x** | 84.6% |
| HCC1395_DORADO | 28,361 | 4,918 | 232 | **21.2x** | 81.8% |
| H2009 | 124,820 | 16,334 | 903 | **18.1x** | 86.2% |
| HCC1937 | 11,826 | 1,146 | 64 | **17.9x** | 89.8% |
| HCC1395 | 28,091 | 3,898 | 242 | **16.1x** | 85.2% |

> **來源**：`q1c_same_locus_concordance_corrected.tsv`

![Enrichment Paired vs TO](figures/fig01_enrichment_paired_vs_to.png)

### 2.3 解讀

**前提**：如果 TO 和 Paired 的 LOH 判定是對等的，`TO_only_LOH` 和 `Paired_only_LOH` 的數量應大致相當。

**證據**：在完全相同的基因組位點上，TO 單獨判定為 LOH 的數量是 Paired 單獨判定為 LOH 的 **16-52 倍**。所有 7 個樣本一致呈現此偏向。

**結論**：TO phasing 系統性地**過度判定 LOH**。這不是隨機波動，而是 TO phasing 機制的固有特性。Concordance 率 81.8-89.8% 看似很高，但不一致的部分中 95%+ 是 TO 多判。

---

## 3. 機制假說：逐步推導

### 3.1 第一步：TO Phasing 機制的根本差異

**前提**：Paired 模式使用 normal sample 提供 germline heterozygous SNV 做 phasing anchor（LongPhase-s），TO 模式沒有 normal sample，只能使用腫瘤自身的 heterozygous SNV 做 phasing anchor（LongPhase-to）。

**關鍵推論**：在 TO 模式下，somatic variant 本身就可能成為 phasing anchor 的一部分。當一個位點被 caller 判定為 somatic variant 並通過 PASS filter，它的 variant allele 天然偏向特定 haplotype。

### 3.2 第二步：Somatic Allele 造成 HP Imbalance

**前提**：真正的 somatic variant (TP) 在腫瘤中具有 clonal allele fraction——支持 variant 的 reads 天然集中在攜帶突變的 haplotype 上。

**證據**：TO TP 的 extreme LOH 比例 = **44.6%**，而 TO FP = 35.9%。TP 的 HP imbalance 顯著高於 FP。

**推導**：TP 是真正的體細胞變異 → variant-supporting reads 偏向一個 haplotype → HP_Ratio 偏離 0.5 → 被判為 LOH-like。因為 TP 的 variant allele 有真實的 haplotype preference，所以 TP 比 FP 更容易呈現 HP imbalance → 更多 TP 被判為 LOH-like → TO LOH 是 TP-enriched。

### 3.3 第三步：TO FP 的 HP Balance 相對較好

**前提**：TO FP 多為 germline variant 或 sequencing artifact，不具有腫瘤特異的 haplotype preference。

**證據**：
- TO FP balanced (HP_Ratio 0.3-0.7) = **41.5%**
- TO TP balanced = **39.7%**
- TO FP 的 HP balance 略好於 TO TP

**推導**：FP 不是真正的 somatic mutation → reads 在 haplotype 間分佈更均勻 → HP_Ratio 更接近 0.5 → 較少被判為 LOH-like。

### 3.4 第四步：Paired FP 的不同機制

**前提**：Paired FP 極少（僅佔 1.04%），且 Paired phasing 由 normal sample 提供，LOH 判定更可靠。

**證據**：Paired FP 的 `effective_hp` 中位數 = 31-38（遠低於 Paired TP 的 59-72）。

**推導**：Paired FP 多發生在低 HP support 區域 → 低 support 下 HP_Ratio 隨機波動更大 → FP 更容易被判為 LOH-like → Paired LOH 是 FP-enriched。這與 TO 的機制完全不同——TO 是 TP 因 somatic allele 被系統性推向 LOH，Paired 是 FP 因低 support 被隨機推向 LOH。

### 3.5 機制總結

| Mode | LOH Enrichment | 根本原因 |
|------|:---:|------|
| **TO 0.805x (TP-enriched)** | TO phasing 在 TP 位點因 somatic allele 造成系統性 HP imbalance → TP 更多被判為 LOH-like |
| **Paired 1.194x (FP-enriched)** | Paired FP 的 HP support 低 → 隨機 LOH 波動更大 → FP 更容易被判為 LOH-like |

---

## 4. HP Balance 分佈數據

### Table 2. 四象限：Extreme LOH 與 Balanced 比例

| Mode | Truth | Extreme LOH (<0.1 or >0.9) | Balanced (0.3-0.7) | 差異方向 |
|------|:---:|:---:|:---:|------|
| **Paired** | TP | 32.3% | 50.1% | |
| **Paired** | FP | **36.5%** | 39.3% | FP 更 LOH-like |
| **TO** | TP | **44.6%** | 39.7% | **TP 更 LOH-like** |
| **TO** | FP | 35.9% | 41.5% | |

> **來源**：`q1e_hp_balance_by_mode_truth.tsv`

![HP Ratio Distribution by Mode x Truth](figures/fig03_hp_ratio_core_hist.png)

**關鍵觀察**：

1. TO TP 的 extreme LOH 比例（44.6%）是四個分組中最高的——高於 TO FP（35.9%）、Paired TP（32.3%）、Paired FP（36.5%）
2. TO TP 的 balanced 比例（39.7%）是四個分組中最低的——TP 在 TO 模式下最不均衡
3. Paired FP 的 extreme LOH 雖然也高於 Paired TP（36.5% vs 32.3%），但差距（+4.2pp）遠小於 TO 的 TP-FP 差距（+8.7pp）

---

## 5. AF 中位數差異：理解 LOH-like 的不同成因

### Table 3. LOH-like Region 的 AF 中位數

| Mode | Truth | LOH Status | AF Median | 解釋 |
|------|:---:|:---:|:---:|------|
| **TO** | TP | LOH-like | **0.580** | 高 AF somatic variant → 強 HP 偏斜 |
| **TO** | FP | LOH-like | **0.969** | Germline homozygous → 極端 HP 但原因不同 |
| **Paired** | TP | LOH-like | **0.938** | 有 normal 校正 → LOH-like 更多是真實 LOH |
| **TO** | TP | non-LOH | 0.433 | 低 AF → HP 較均衡 |
| **TO** | FP | non-LOH | 0.880 | |

> **來源**：`q1a_af_by_mode_truth_loh.tsv`

![AF Violin by Truth x LOH Status](figures/fig02_af_violin_by_truth_loh.png)

![AF vs HP Ratio Scatter](figures/fig08_af_vs_hp_ratio_scatter.png)

**關鍵推導**：

1. **TO TP LOH-like (AF=0.580)**：中等偏高的 somatic AF → 足夠的 variant-supporting reads 偏向一個 haplotype → HP imbalance → 被判為 LOH-like。**這是 somatic clonality 的自然結果**。
2. **TO FP LOH-like (AF=0.969)**：接近 homozygous 的 germline variant → 幾乎所有 reads 都攜帶 variant → 極端 HP_Ratio。但這些 FP 的「LOH-like」與 somatic LOH 機制不同——它們是 germline homozygosity。
3. **Paired TP LOH-like (AF=0.938)**：有 normal sample 校正後，LOH-like 的 TP 多數 AF 極高 → 更可能是真正的 LOH block 中的體細胞變異。

---

## 6. TO LOH 作為 TP Rescue 候選

### 6.1 低 AF TP 的 LOH-like 比例

既然 TO LOH 是 TP-enriched，一個自然的問題是：能否利用 LOH-like 特徵來 rescue 低可信度的 TP？

#### Table 4. 低 AF TP 的 LOH-like 比例

| AF 閾值 | TO TP 數量 | TO LOH-like% | Paired TP LOH-like% | TO 優勢 |
|:---:|---:|:---:|:---:|:---:|
| < 0.10 | 1,457 | **51.1%** | 45.4% | +5.7pp |
| < 0.15 | 9,435 | **44.3%** | 15.1% | +29.2pp |
| < 0.20 | 21,284 | **40.6%** | 7.1% | +33.5pp |
| < 0.30 | 59,109 | **35.2%** | 8.4% | +26.8pp |

> **來源**：`q2c_low_af_tp_loh_fraction.tsv`

![Low AF TP LOH Rescue Potential](figures/fig09_low_af_tp_loh_rescue.png)

**發現**：AF < 0.10 的 TO TP 有超過一半（51.1%）為 LOH-like。即使在最低 AF 的 TP 中，LOH-like 仍然是強訊號。

**但必須注意的限制**：

1. 這些 TP 已在 caller 的最終 PASS VCF 中（不是 FN）
2. 真正的 rescue 價值在於：**caller 已過濾但實際為 TP 的 FN 位點**是否也有高 LOH-like
3. 目前沒有 FN regions 的 LOH 數據（FN 不在 ISM pipeline 中）→ 需要新的數據收集

### 6.2 LOH-like by VerificationClass

ISM 的 VerificationClass 是基於甲基化距離矩陣的 PERMANOVA 結果分類。如果 LOH-like 的 TP enrichment 在所有 VerificationClass 中都一致，說明這是普遍現象而非特定 class 驅動。

#### Table 5. LOH-like by VerificationClass（TO 模式）

| VerificationClass | TP LOH-like% | FP LOH-like% | TP 優勢 (pp) |
|:---:|:---:|:---:|:---:|
| **Strong** | 33.5% | 19.8% | **+13.7pp** |
| **Weak** | 30.8% | 17.3% | **+13.5pp** |
| **Noise** | 56.1% | 50.5% | +5.6pp |
| **Subclone** | 59.4% | 55.4% | +4.0pp |

> **來源**：`q2a_loh_by_verification_class.tsv`

**結論**：所有 VerificationClass 中 TO TP 的 LOH-like 比例都高於 FP，無例外。Strong 和 Weak class 的差距最大（+13.5-13.7pp），Noise 和 Subclone 差距較小但方向一致。

### 6.3 TO TP LOH vs non-LOH 品質差異

#### Table 6. TO TP LOH-like vs non-LOH 品質指標

| 指標 | TO TP LOH-like | TO TP non-LOH | 差異 |
|------|:---:|:---:|------|
| AF 中位數 | 0.580 | 0.433 | LOH TP 的 AF 更高 |
| Quality Score 中位數 | 75 | 75 | 無差異 |
| Effective HP 中位數 | 62 | 79 | LOH TP 的 HP support 較低 |
| HP0 ratio 中位數 | 0.0 | 0.0 | 無差異 |
| Tier A/A+ 比例 | 87-89% | 89-93% | 接近 |

> **來源**：`q2b_tp_loh_vs_noloh_quality.tsv`

LOH-like TP 和 non-LOH TP 在品質指標上差異不大，最主要差異在 AF（0.580 vs 0.433）。LOH-like TP 的 AF 顯著更高，進一步支持「somatic allele 的 clonal fraction 驅動 HP imbalance」的機制假說。

---

## 7. 本章核心結論

### 確認的結論

| # | 結論 | 強度 | 關鍵數據 |
|---|------|:---:|------|
| 1 | TO phasing 系統性過判 LOH：同位點 TO_only_LOH 是 paired 的 16-52x | 確定 | `q1c_same_locus_concordance_corrected.tsv`，288,609 loci |
| 2 | 過判原因：TO phasing 中 somatic allele 造成系統性 HP imbalance | 強支持 | HP balance 四象限 + AF 中位數差異 |
| 3 | TO LOH-like 在所有 VerificationClass 中都是 TP > FP | 確定 | `q2a_loh_by_verification_class.tsv` |
| 4 | 低 AF TP (< 0.1) 有 51.1% 為 LOH-like，具有 rescue 潛力 | 初步 | `q2c_low_af_tp_loh_fraction.tsv`（但缺 FN 驗證） |
| 5 | Paired FP LOH-enriched 原因：FP effective_hp 低 → 隨機波動 → 被判 LOH-like | 強支持 | effective_hp 中位數差異（FP 31-38 vs TP 59-72） |

### 不能宣稱的

1. **LOH-like 可以直接用於 TP rescue**——缺少 FN 數據驗證。目前只能說 PASS VCF 中的低 AF TP 有高 LOH-like，但真正的 rescue 需要 FN 位點的 LOH 數據。
2. **TO phasing 的 LOH-like 都是「真實 LOH」**——很大一部分是 phasing 層面的 HP imbalance（somatic allele 造成），不一定是基因組層面的 LOH event。
3. **機制假說已完全驗證**——第 3 節的推導是基於間接證據（AF、HP balance、concordance），尚未有直接的 per-read haplotype assignment 驗證。

---

## 8. 資料來源與可追溯性

| 項目 | 路徑 |
|------|------|
| 分析腳本 | `scripts/analysis/build_post_hp_fix_to_loh_investigation.py` |
| Workspace | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_post_hp_fix_to_loh_investigation/` |
| 同位點 concordance | `q1c_same_locus_concordance_corrected.tsv` |
| AF by mode/truth/LOH | `q1a_af_by_mode_truth_loh.tsv` |
| HP balance 四象限 | `q1e_hp_balance_by_mode_truth.tsv` |
| VerificationClass LOH | `q2a_loh_by_verification_class.tsv` |
| TP LOH vs non-LOH 品質 | `q2b_tp_loh_vs_noloh_quality.tsv` |
| 低 AF TP LOH fraction | `q2c_low_af_tp_loh_fraction.tsv` |
| Master dataset | `20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz` |
| Pre-fix dataset | `20260327_loh_round1_cross_sample_audit_before_hp_fix/all_region_rows.tsv.gz` |
| HP bug 修正 | `src/core/ReadParser.cpp` switch-case mapping（HP:i:11->"1-1", HP:i:21->"2-1", HP:i:33->"3"） |

---

## 待驗證問題（已驗證 / 已更新）

### 已解決

4. **~~LOH-like TP 的 AF=0.580 是否是 clonal fraction proxy~~** ✅ **是的。** Q3 機制分析確認：TP LOH rate +51.6% 由 somatic allele 的 haplotype preference 驅動，AF=0.580 反映 clonal fraction。TO 的「LOH-like」本質上偵測的是 somatic clonality 而非基因組 LOH。（見 Section 02 C13）

5. **~~COLO829 concordance 25.1x 是否受低 eff_hp 影響~~** ✅ **是的。** Q1 根因分析確認：COLO829 低測序深度 (~30x, median num_reads=29) 導致 hp_ratio 離散化嚴重，concordance 比值被放大。（見 Section 03 C11）

### 部分解決

2. **區分「真實 LOH」vs「phasing HP imbalance」** — Q3 分析提供量化：39,724 個從 Paired non-LOH 翻轉為 TO LOH 的位點中，**71.6% 在 TO 的 min(HP1,HP2)=0**（完全單 haplotype），同位點在 Paired 維持 HP1:HP2=8:8。GT completeness / PS block size 的進一步精細分離仍需新分析。

### 尚未解決

1. **FN 位點的 LOH-like 比例**（P0 優先）：需 FN ISM 數據。對應 Section 09 P2「O9 FN 觀察」。
3. **Per-read haplotype assignment 直接驗證**：需逐 read 層級檢查，超出現有數據範圍。定位 P2。

---

## 認知門檻補充建議

1. **TO phasing 與 Paired phasing 的根本差異**：Paired 模式使用 normal sample 的 germline heterozygous SNV 做 phasing anchor（LongPhase-s），提供穩定的 haplotype scaffold。TO 模式只能使用腫瘤自身的 heterozygous variant（可能包含 somatic variant）做 anchor（LongPhase-to）。當 somatic variant 本身成為 phasing anchor 時，variant-supporting reads 天然集中在一個 haplotype，導致 HP_Ratio 系統性偏離 0.5。

2. **「LOH-like」不等於「基因組 LOH」**：在 ISM 中，「LOH-like」的操作定義是 `HP_Ratio < 0.1 or > 0.9`（即 extreme haplotype imbalance）。在 TO 模式下，這個 imbalance 的主要成因是 somatic allele 的 haplotype preference（phasing artifact），而非基因組層面的 LOH event（一個 haplotype 的實際缺失）。因此 TO 的「LOH-like」更準確的名稱應該是「HP-imbalanced」。

3. **Enrichment 方向翻轉的實務意義**：在 Paired 模式下，LOH 區域是 FP 高風險區（LOH penalty 有道理）。在 TO 模式下，LOH 區域反而是 TP 富集區（LOH penalty 反向作用）。這意味著 **Paired 和 TO 不能使用同一套 Quality Score 公式**——至少 LOH penalty 項必須 mode-specific。

4. **同位點 concordance 的統計意義**：16-52x 的比值不代表 TO 判錯了 16-52 倍。85% concordance 代表大多數位點判定一致。比值的意義在於：**不一致的部分幾乎全是 TO 多判**（而非 Paired 多判），確認了系統性偏向的方向。

5. **TP rescue 的前提條件**：利用 LOH-like 做 TP rescue 的邏輯是「LOH-like 在 TP 中更常見 → 如果一個被 caller 過濾的位點是 LOH-like，它更可能是 TP」。但這需要兩個前提：(a) FN 位點確實有高 LOH-like 比例，(b) LOH-like 的 TP 優勢在 FN 群體中仍然成立。兩者目前都缺乏驗證數據。
