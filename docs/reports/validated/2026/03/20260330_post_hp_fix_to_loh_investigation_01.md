<!--
建立時間: 2026-03-30 21:00
目標: HP integer tag bug 修正後 TO LOH 深度調查：enrichment 方向差異原因、TP rescue 潛力、Tier 結構重塑
處理範圍: Q1 為什麼 TO 和 paired LOH enrichment 方向相反、Q2 TO LOH TP rescue、Q3 Tier 結構視覺化
關聯檔案:
  - scripts/analysis/build_post_hp_fix_to_loh_investigation.py
  - docs/reports/validated/2026/03/20260330_TO_LOH_enrichment_post_hp_fix_01.md
  - docs/architecture/20260327_InterSubMod研究願景定錨_01.md
  - docs/plans/2026/03/20260327_LOH_round2_execution_spec_01.md
workspace: /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_post_hp_fix_to_loh_investigation
-->

# Post-HP-Fix TO LOH 深度調查報告

> **背景**：2026-03-30 修正 ReadParser.cpp HP integer tag bug 後，TO LOH enrichment 確認為 **0.805**（TP 富集），與 paired 的 **1.194**（FP 富集）方向相反。本報告深入調查原因、評估 rescue 潛力、並視覺化 Tier 結構重塑。

---

## 1. 核心發現摘要

| 發現 | 數據依據 | 影響 |
|------|---------|------|
| **TO phasing 在同位點產生的 LOH-like 是 paired 的 16-52 倍** | `q1c_same_locus_concordance_corrected.tsv` | 這是 enrichment 方向差異的根本原因 |
| **TO TP 極端 LOH 比例 44.6% >> TO FP 35.9%** | `q1e_hp_balance_by_mode_truth.tsv` | TO LOH 是 TP marker，不適合作 FP filter |
| **低 AF TP（<0.1）有 51.1% 為 LOH-like** | `q2c_low_af_tp_loh_fraction.tsv` | LOH-like 可作為低 AF TP 的 rescue 信號 |
| **修正後 88% TO regions 在 Tier A/A+** | `q3a_tier_distribution_before_after.tsv` | TO phasing 品質遠比之前認知的好 |

---

## 2. Q1：為什麼 TO 和 paired LOH enrichment 方向相反？

### 2.1 同位點 LOH concordance（最關鍵證據）

> 資料：`q1c_same_locus_concordance_corrected.tsv`，288,609 個 valid 同位點 loci
> 圖表：`fig01_enrichment_paired_vs_to.png`

| 樣本 | 同位點數 | TO_only_LOH | paired_only_LOH | 比值 | concordance |
|------|---------|------------|-----------------|------|------------|
| HCC1395 | 28,091 TP | 3,898 | 242 | **16.1×** | 85.2% |
| HCC1395_DORADO | 28,361 TP | 4,918 | 232 | **21.2×** | 81.8% |
| COLO829 | 32,809 TP | 5,205 | 207 | **25.1×** | 83.5% |
| H1437 | 44,843 TP | 5,802 | 112 | **51.8×** | 86.8% |
| H2009 | 124,820 TP | 16,334 | 903 | **18.1×** | 86.2% |
| HCC1937 | 11,826 TP | 1,146 | 64 | **17.9×** | 89.8% |
| HCC1954 | 16,342 TP | 2,421 | 96 | **25.2×** | 84.6% |

**結論**：在完全相同的基因組位點上，TO phasing（longphase-to）產生的 LOH-like 判定比 paired phasing（longphase-s）多 16-52 倍。這不是隨機波動，而是系統性差異。

### 2.2 機制解釋

> 資料：`q1a_af_by_mode_truth_loh.tsv`、`q1e_hp_balance_by_mode_truth.tsv`
> 圖表：`fig02_af_violin_by_truth_loh.png`、`fig03_hp_ratio_core_hist.png`、`fig08_af_vs_hp_ratio_scatter.png`

**HP balance 分佈（關鍵數據）**：

| Mode | Truth | Extreme LOH (< 0.1 or > 0.9) | Balanced (0.3-0.7) | 差異方向 |
|------|-------|------|---------|---------|
| **Paired** | TP | 32.3% | 50.1% | FP 更 LOH-like |
| **Paired** | FP | **36.5%** | 39.3% | ↑ |
| **TO** | TP | **44.6%** | 39.7% | **TP 更 LOH-like** |
| **TO** | FP | 35.9% | 41.5% | ↓ |

**解釋假說**（基於數據支持）：

1. **TO phasing 在 TP 位點的 somatic allele 造成 HP imbalance**
   - TO TP LOH-like AF 中位數 = 0.580（高 AF somatic variant → 強 HP 偏斜）
   - TO FP LOH-like AF 中位數 = 0.969（germline homozygous → HP 極端但原因不同）
   - Paired TP LOH-like AF 中位數 = 0.938（paired 有 normal 校正，LOH-like 更多是真實 LOH）

2. **TO FP 多為 germline variant，HP balance 相對好**
   - TO FP balanced (0.3-0.7) = 41.5% > TO TP balanced = 39.7%
   - Paired FP balanced = 39.3%（paired FP 本身就少，且 HP support 低 → 更容易 LOH-like 波動）

3. **Paired FP 的 HP support 低 → 隨機 LOH fluctuation**
   - Paired FP effective_hp 中位數 = 31-38（遠低於 TP 的 59-72）
   - 低 support + 隨機波動 → FP 更容易被判為 LOH-like → FP enrichment

### 2.3 小結

| Mode | LOH enrichment | 主要原因 |
|------|---------------|---------|
| **Paired 1.194× (FP enriched)** | FP 的 HP support 低 → 隨機 LOH 波動更大 → FP 更容易被判 LOH-like |
| **TO 0.805× (TP enriched)** | TO phasing 在 TP 位點的 somatic allele 造成系統性 HP imbalance → TP 更多 LOH-like |

---

## 3. Q2：TO LOH 作為 TP rescue 候選

### 3.1 低 AF TP 的 LOH-like 比例

> 資料：`q2c_low_af_tp_loh_fraction.tsv`
> 圖表：`fig09_low_af_tp_loh_rescue.png`

| AF 閾值 | TO TP 數量 | LOH-like 比例 | Paired TP LOH-like 比例 |
|---------|-----------|--------------|----------------------|
| < 0.10 | 1,457 | **51.1%** | 45.4% |
| < 0.15 | 9,435 | **44.3%** | 15.1% |
| < 0.20 | 21,284 | **40.6%** | 7.1% |
| < 0.30 | 59,109 | **35.2%** | 8.4% |

**解讀**：低 AF 的 TO TP 有極高比例為 LOH-like（AF < 0.1 時超過一半）。這暗示 LOH-like 特徵可以作為低可信度 TP 的補充 evidence，但需注意：

- 這些 TP 已在 caller 的最終 PASS VCF 中（不是 FN）
- 真正的 rescue 價值在於：**caller 已過濾但實際為 TP 的 FN 位點**是否也有高 LOH-like
- 目前沒有 FN regions 的 LOH 數據（FN 不在 ISM pipeline 中）→ 需要新的數據收集

### 3.2 LOH-like by VerificationClass

> 資料：`q2a_loh_by_verification_class.tsv`

| VerificationClass | TP LOH-like | FP LOH-like | TP 更高？ |
|-------------------|-------------|-------------|----------|
| **Noise** | 56.1% | 50.5% | ✅ +5.6pp |
| **Strong** | 33.5% | 19.8% | ✅ +13.7pp |
| **Subclone** | 59.4% | 55.4% | ✅ +4.0pp |
| **Weak** | 30.8% | 17.3% | ✅ +13.5pp |

所有 VerificationClass 中 TO TP 的 LOH-like 都高於 FP，尤其在 Strong 和 Weak class 差距達 13+ 百分點。

### 3.3 TO TP LOH vs non-LOH 品質差異

> 資料：`q2b_tp_loh_vs_noloh_quality.tsv`
> 圖表：`fig07_to_tp_loh_vs_noloh_features.png`

| 指標 | TO TP LOH-like | TO TP non-LOH | 差異 |
|------|---------------|---------------|------|
| AF 中位數 | 0.580 | 0.433 | LOH TP 的 AF 更高 |
| Quality Score 中位數 | 75 | 75 | 無差異 |
| Effective HP 中位數 | 62 | 79 | LOH TP 的 HP support 較低 |
| HP0 ratio 中位數 | 0.0 | 0.0 | 無差異 |
| Tier A/A+ 比例 | 87-89% | 89-93% | 接近 |

**解讀**：LOH-like TP 和 non-LOH TP 在品質指標上差異不大，但 LOH-like TP 的 AF 顯著更高（0.580 vs 0.433），暗示 LOH-like 與 clonal allele fraction 有關。

---

## 4. Q3：HP Fix 後 Tier 結構重塑

### 4.1 TO Tier 分佈

> 資料：`q3a_tier_distribution_before_after.tsv`
> 圖表：`fig04_tier_before_after.png`、`fig06_per_sample_tier_shift.png`

| Tier | Before Fix | After Fix | 變化 |
|------|-----------|-----------|------|
| **A+** (≥50) | 29.5% (123,619) | **71.8%** (301,151) | +143% |
| **A** (30-49) | 21.4% (89,807) | 16.5% (69,432) | -23% |
| B (10-29) | 24.4% (102,461) | 10.0% (41,928) | -59% |
| C (1-9) | 15.3% (64,006) | 1.5% (6,111) | -90% |
| C0 (0) | 9.5% (39,799) | **0.3%** (1,070) | -97% |

**認知修正**：修正前以為 TO 有 ~25% regions 在低品質 tier（C/C0），實際上只有 ~2%。TO phasing 品質遠好於之前的認知。

### 4.2 Enrichment by Tier 變化

> 資料：`q3c_enrichment_by_tier_before_after.tsv`
> 圖表：`fig05_enrichment_by_tier_before_after.png`

| Mode | Tier | Before Fix Enrichment | After Fix Enrichment | 變化 |
|------|------|---------------------|---------------------|------|
| TO | A+ | 0.815 | **0.770** | TP 富集更強 |
| TO | A | 0.962 | 0.773 | 從接近中性變為 TP 富集 |
| TO | B | 1.002 | **0.938** | 從中性變為弱 TP 富集 |
| TO | C | 0.990 | 1.015 | 接近中性（樣本少） |
| Paired | A+ | 1.236 | 1.236 | 不受影響 |
| Paired | A | 0.948 | 0.948 | 不受影響 |

**修正前的 TO Tier B/C 包含了大量被 bug 低估的 regions**，修正後這些 regions 升到 A/A+，使得：
- A/A+ tier 的 TP 富集更強（更多正確分配的 TP regions）
- B/C tier 變得非常小（剩下的才是真正低品質的 regions）

### 4.3 HP0 by Tier

> 圖表：`fig10_hp0_by_tier_violin.png`

修正後 HP0 ratio 呈清晰的單調遞減：
- C0: ~85% HP0（幾乎全是 unphased）
- C: ~78% HP0
- B: ~22% HP0
- A: ~9% HP0
- A+: ~3% HP0

這個模式在修正前是混亂的（因為大量 regions 被錯誤歸類），現在非常合邏輯。

---

## 5. 圖表索引

| 圖表 | 描述 | 對應問題 |
|------|------|---------|
| `fig01_enrichment_paired_vs_to.png` | 7 樣本 enrichment 雙向對比 | Q1 |
| `fig02_af_violin_by_truth_loh.png` | AF 分佈 by truth × LOH status | Q1 |
| `fig03_hp_ratio_core_hist.png` | HP ratio 直方圖 by mode × truth | Q1 |
| `fig04_tier_before_after.png` | Tier 結構 before/after 堆疊圖 | Q3 |
| `fig05_enrichment_by_tier_before_after.png` | Enrichment by tier 折線圖 | Q3 |
| `fig06_per_sample_tier_shift.png` | 7 樣本 tier 變化 | Q3 |
| `fig07_to_tp_loh_vs_noloh_features.png` | TO TP 特徵比較 | Q2 |
| `fig08_af_vs_hp_ratio_scatter.png` | AF vs HP ratio 散布圖 | Q1 |
| `fig09_low_af_tp_loh_rescue.png` | 低 AF TP 的 LOH rescue 潛力 | Q2 |
| `fig10_hp0_by_tier_violin.png` | HP0 ratio by tier violin | Q3 |

---

## 6. 結論與下一步

### 確認的結論

1. **TO LOH enrichment 與 paired 方向相反的根本原因**：TO phasing 在 TP 位點因 somatic allele 產生系統性 HP imbalance，同位點 TO_only_LOH 是 paired 的 16-52 倍
2. **TO LOH-like 是 TP marker**：所有 VerificationClass 中 TO TP 的 LOH-like 比例都高於 FP
3. **TO phasing 品質遠優於之前認知**：88% regions 在 Tier A/A+，只有 2% 在 C/C0
4. **低 AF TP 有高比例 LOH-like**：AF < 0.1 的 TO TP 有 51.1% 為 LOH-like

### 不能宣稱的

- LOH-like 可以直接用於 TP rescue（缺少 FN 數據驗證）
- TO phasing 的 LOH-like 都是「真實 LOH」（很大一部分是 phasing 層面的 HP imbalance，不一定是基因組層面的 LOH）

### 建議下一步

| 優先度 | 方向 | 需要的資料 |
|--------|------|-----------|
| **P0** | 取得 FN regions 的 LOH 數據，驗證 LOH rescue 可行性 | 需要從 benchmark_comparison.tsv 中提取 FN 位點，跑 ISM pipeline |
| **P1** | 區分「真實 LOH」vs「phasing HP imbalance」| 同位點 LOH subtype 分析、PS block boundary 效應 |
| **P2** | 將 LOH TP enrichment 整合進 Phase 1 ML feature | 將 `core_loh_like` 或 `hp_ratio_core` 作為 read classifier 輸入特徵 |

---

## 7. 資料來源與可追溯性

| 項目 | 路徑 |
|------|------|
| 分析腳本 | `scripts/analysis/build_post_hp_fix_to_loh_investigation.py` |
| Workspace | `output/synthesis/observation_workspaces/20260330_post_hp_fix_to_loh_investigation/` |
| Master dataset | `output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz` |
| Pre-fix dataset | `output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit_before_hp_fix/all_region_rows.tsv.gz` |
| HP bug 修正 commit | ReadParser.cpp switch-case mapping（HP:i:11→"1-1", HP:i:21→"2-1", HP:i:33→"3"） |
