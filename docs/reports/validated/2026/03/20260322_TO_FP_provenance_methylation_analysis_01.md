<!--
建立時間: 2026-03-22
目標: HCC1395 5kHz TO FP 五層分類 + 甲基化特徵分析 + F1 提升可行性評估
處理範圍: HCC1395 5kHz TO track 全量 FP (11,598) + ISM 候選池 (298 FP + 773 TP)
關聯檔案:
  - research/fp_provenance/20260322_hcc1395_5khz_to/
  - scripts/analysis/build_to_fp_provenance_analysis.py
  - /big8_disk/.../20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_joined_features.tsv
  - /big8_disk/.../20260311_phase2_annotation_layer/annotated_candidates/hcc1395_5khz_to_annotated_candidates.tsv
-->

# TO FP Provenance 分類與甲基化特徵分析

> 生成時間：2026-03-22
> 資料集：HCC1395 5kHz TO (5kHz_simplex_5mCG_5hmCG)
> 前置研究：autoresearch v1.1.0（H012 GQ>=3 已達最佳 rescue F1 +0.009365）

---

## 一、研究背景與動機

InterSubMod TO track 目前最佳 F1：
- 基線 F1 = 0.712697（LongPhase-TO 後）
- H012 rescue 後 F1 = ~0.722062（GQ>=3 rescue +0.009365）

TO 有 11,843 個 FP（~18.9x 多於 paired 627 個 FP）。本研究試圖回答：

**研究問題**：在 TO 的 FP 中，哪些可以透過 normal 資料解決、哪些是 TO 特有？ISM 甲基化特徵是否能進一步區分 TP 與 FP，以提升 TO F1？

---

## 二、FP 五層過濾分類結果

### 2.1 全量 FP 統計（HCC1395 5kHz TO）

| 類別 | 數量 | 佔比 | 說明 |
|------|------|------|------|
| **caller_pon_filtered** | 2,717,339 | 99.5% | NonSomatic/PON 過濾（正確） |
| **caller_nonpon_filtered** | 2,596 | 0.095% | LowQual/VariantCluster 等過濾（正確） |
| **longphase_to_removed** | 0 | 0% | LongPhase-TO 無獨立 FP 過濾 |
| **to_postprocess_removed** | 8 | <0.001% | ISM 後處理規則移除（AF<0.15 & AD>0.15） |
| **to_residual_final_fp** | **11,598** | 0.425% | **進入最終輸出的錯誤 calls（核心研究對象）** |
| 合計 | 2,731,541 | 100% | 所有在 truth region 內的 FP |

**關鍵發現 1**：絕大多數 FP（99.5%）已被 PON/NonSomatic 正確過濾，且均為高品質過濾（gnomAD/dbSNP/1000G/CoLoRSdb 命中）。

**關鍵發現 2**：LongPhase-TO 的 NoAncestry/MultiHap flags 在此資料集上**對 FP 無獨立過濾貢獻**（longphase_to_removed = 0）。

### 2.2 殘餘 FP（11,598）的 Paired Oracle 解析

| Paired 解析方式 | 數量 | 佔 11,598% | 說明 |
|----------------|------|-----------|------|
| **raw_absent** | 11,430 | **98.6%** | Paired VCF 根本不輸出此 variant（最大群） |
| **persistent** | 87 | 0.75% | 即使有 normal 也無法解決（硬性 FP） |
| **raw_filter** | 63 | 0.54% | Paired VCF 中被 FILTER 標記 |
| **longphase_s** | 18 | 0.15% | Paired LongPhase-S phasing 移除 |

**最重要發現**：**98.6% 的殘餘 FP（11,430 個）在 Paired 模式下根本不被 ClairS 輸出**。這意味著這些 FP 是 ClairS-TO 特有的 calling 行為所產生，不是簡單的「過濾」問題，而是「calling model 差異」問題。

**僅有 87 個 FP（0.75%）是即使有 normal 也無法解決的持久性 FP**。

---

## 三、ISM 候選池 FP 分類分析

### 3.1 ISM 298 FP 候選的歸屬

| ISM FP 候選歸屬 | 數量 | 佔 298% | 說明 |
|----------------|------|---------|------|
| **to_residual_final_fp** | 208 | 69.8% | 仍在最終 PASS 輸出中（目標過濾對象） |
| **caller_nonpon_filtered** | 85 | 28.5% | 已被 caller 正確過濾（不在輸出中） |
| **to_postprocess_removed** | 5 | 1.7% | 已被後處理移除（不在輸出中） |

### 3.2 208 個 PASS FP 候選的 Paired 解析

| Paired 解析方式 | 數量 | 佔 208% |
|----------------|------|---------|
| raw_absent（TO 特有問題）| 148 | 71.2% |
| persistent（硬性 FP）| 26 | 12.5% |
| raw_filter（Paired 可過濾）| 23 | 11.1% |
| longphase_s（Phasing 可移除）| 11 | 5.3% |

**結論**：ISM 候選池中，208 個 PASS FP 裡，**174 個（83.7%）是 TO-hard FP**（raw_absent + persistent），只有 34 個（16.3%）在有 normal 資料時可解決。

---

## 四、甲基化特徵差異分析

### 4.1 ISM 候選池：TP vs 各 FP 群體特徵比較

| 特徵 | TP (n=773) | FP_raw_absent (n=210) | FP_persistent (n=31) | FP_raw_filter (n=44) |
|------|-----------|----------------------|---------------------|----------------------|
| Quality_Score (中位數) | **75.0** | 60.0 | 75.0 | 50.0 |
| PairwiseMedianDist | **0.210** | 0.175 | 0.217 | 0.114 |
| CramersV | 0.000 | 0.000 | 0.000 | 0.000 |
| AlleleDelta | 0.010 | 0.000 | 0.000 | 0.000 |
| hp_assign_rate | **0.963** | 0.928 | 0.925 | 0.420 |
| GQ (中位數) | **14** | 8 | 10 | 4 |
| AF (中位數) | **0.214** | 0.241 | 0.147 | 0.105 |

### 4.2 Mann-Whitney U 統計顯著性

| 比較群體 | QS | PMD | GQ | AF | CV |
|---------|----|----|----|----|-----|
| TP vs FP_raw_absent | *** | ** | *** | * | ns |
| TP vs FP_persistent | ns | ns | * | *** | ns |
| TP vs FP_raw_filter | *** | *** | *** | *** | ns |
| TP vs FP_caller_nonpon | *** | *** | *** | *** | * |

**重要觀察**：
- **FP_raw_absent**：GQ 和 QS 有顯著差異（GQ:8 vs TP:14，QS:60 vs 75）
- **FP_persistent**：AF 是唯一顯著特徵（AF:0.147 vs TP:0.214，低 VAF FP）
- **CramersV** 在所有 FP 群體中均無顯著差異（TO HP-based 特徵根本限制）

---

## 五、全量 PASS 輸出的特徵判別力分析

### 5.1 GQ 分佈（全量 PASS output）

| 百分位數 | TP (n=28,509) | FP (n=11,598) |
|---------|--------------|--------------|
| P10 | 14 | 15 |
| P25 | 17 | 17 |
| P50 | 19 | **19** |
| P75 | 21 | 21 |
| P90 | 24 | 24 |

**GQ 分佈幾乎完全相同** → GQ 過濾 PASS 輸出無效。

GQ < 10 filter 掃描結果：FP/TP 移除比 = 0.34（每移除 1 個 FP，損失約 2.9 個 TP）→ **delta_F1 = -0.008628**（負增益）。

### 5.2 QUAL 分佈（全量 PASS output）

| 百分位數 | TP | FP |
|---------|----|----|
| P50 | 19.6 | 19.9 |
| P75 | 21.6 | 21.9 |

**QUAL 分佈幾乎完全相同** → QUAL 過濾無效。

### 5.3 Quality_Score 分佈（全量 ISM 分析後 PASS output）

| 閾值 | TP (n=28,495) | FP (n=11,601) | 比例 |
|-----|--------------|--------------|------|
| QS < 50 | 7.4% | 7.4% | 完全一致 |
| QS < 60 | 22.2% | 22.1% | 幾乎完全一致 |
| QS < 75 | 34.2% | 34.7% | 幾乎完全一致 |

**Quality_Score 在全量 PASS 輸出中完全無法區分 TP 與 FP**（各百分位數差異 < 0.5%）。

### 5.4 AF 分佈分析

AF 高端（0.7-0.9）確實有 FP 富集：
- AF 0.8-0.9：FP fraction = **52%**（基線整體 FP fraction = 28.9%）
- AF 0.7-0.8：FP fraction = **45.5%**

但由於高 AF 範圍同時包含大量 LOH 型態的 TP（高純度腫瘤體細胞 SNV），AF 過濾仍為負增益：
- AF > 0.7 filter → FP/TP 移除比 = 0.88（每移除 1 FP，損失 1.14 TP）→ **delta_F1 = -0.057**

---

## 六、核心結論

### 6.1 TO FP 問題的根本性質

```
TO FP 問題的根源：「ClairS-TO 在無 normal 下的 calling 行為差異」

11,430/11,598 (98.6%) 殘餘 FP 屬於 "raw_absent" 類別：
→ 這些 FP 在 Paired 模式下根本不被 ClairS 輸出（連 FILTER 都不會出現）
→ 這不是「過濾問題」，而是「calling model 差異」

意義：
(a) 有 normal 資料 = ClairS 自然不會 call 這些 variants（信號被 normal 抑制）
(b) 無 normal (TO) = ClairS 看到「疑似 somatic 信號」但無法確認
(c) ISM 甲基化特徵無法替代 normal 樣本的信息量
```

### 6.2 甲基化特徵判別力結論

| 分析範圍 | ISM 特徵判別力 | 結論 |
|---------|--------------|------|
| 全量 PASS output（28,509 TP + 11,598 FP）| **完全無判別力** | QS/GQ/QUAL 分佈完全相同 |
| ISM 低 GQ 候選池（773 TP + 298 FP）| **部分判別力** | QS 和 GQ 在低品質範圍有差異 |
| Persistent FP（87 個）| AF 有弱信號 | AF 低於 TP（0.147 vs 0.214），但 n=87 太小 |

### 6.3 F1 提升可行性評估

| 方法 | 預期 delta_F1 | 可行性 |
|------|--------------|-------|
| GQ < threshold 過濾（任意閾值）| **負值** | ❌ 不可行 |
| QS < threshold 過濾（任意閾值）| **負值** | ❌ 不可行（TP/FP 分佈相同）|
| AF 區間過濾 | **負值** | ❌ 不可行 |
| Persistent FP AF < 0.15 過濾 | +0.001~0.002 | ⚠️ 極小（僅 87 個 FP）|
| 擴大 ISM 覆蓋率（7%→50%）| 未知但潛力高 | ✅ 最有潛力（FN rescue 方向）|
| H012 (GQ>=3) rescue（已採納）| +0.009365 | ✅ 最有效 |

---

## 七、研究局限與重要警告

1. **ISM FP 候選池代表性問題**：
   - ISM 分析的 298 FP candidates 僅是全量 11,598 FP 的 2.5%
   - 且是 **低 GQ 偏向樣本**（ISM 分析範圍優先低品質邊界 variants）
   - 不能由此外推全量 FP 的 ISM 特徵行為

2. **TO haplotagging 根本限制**（已在 H001-H012 autoresearch 確認）：
   - AlleleDelta/HPP/AlleleMethDelta 在 TO 全數無效
   - CramersV = 0 for ALL groups → HP-based 特徵不可用

3. **覆蓋率缺口**：
   - ISM 目前只分析了 7% 的 FN pool（773/11,051）
   - 擴大 ISM 覆蓋率是提升 rescue 潛力的關鍵，非 FP 過濾

4. **Persistent FP（87 個）**：
   - AF 信號存在但樣本量過小（87 個）
   - 即使完美過濾全部 87 個（不損失任何 TP）：delta_F1 = +0.0015
   - 實際效益受限，不建議作為優先方向

---

## 八、下一步建議

### P1（最高優先）：擴大 ISM FN 覆蓋率（rescue 方向）
```
當前：ISM 覆蓋 773/11,051 FN = 7%
目標：提升到 30-50%（覆蓋 3,315-5,526 FN）
潛力：若 rescue precision ~74%（H012 水準），可獲得 ~2,453-4,090 額外 TP
預估 delta_F1：+0.03 至 +0.06（理論最大值範圍）

方法：
- 降低 ISM region selection 門檻
- 對更大範圍的 FN candidates 執行 ISM 分析
```

### P2（次要）：DORADO TO 上驗證 H012
```
已知 DORADO TO rescue 空間極小（+0.000476 max, Phase 2）
目的：確認 GQ>=3 不在 DORADO TO 引入大量 FP
方法：執行同樣的 GQ sweep 在 DORADO TO 資料
```

### P3（觀察性）：Persistent FP AF 分析
```
87 個 persistent FP（即使 paired 也無法解決）
AF 中位數 = 0.147 vs TP 0.214（p<0.001）
建議：確認這些 FP 的生物學特性（是否為 MNP/複雜變異？）
注意：過濾 87 個 FP 的最大收益 delta_F1 ≈ +0.0015（不建議優先投資）
```

---

## 九、方法論貢獻

本研究建立了 TO FP 五層分類框架：

```
層級  過濾機制              FP 數量    特性
─────────────────────────────────────────────────────────────
L1   PON/NonSomatic 過濾   2,717,339   99.5% 有效，ClairS 已處理
L2   非 PON Caller 過濾    2,596        0.095%，LowQual/VariantCluster
L3   LongPhase-TO 過濾     0            本資料集無貢獻
L4   ISM 後處理規則         8            極少，AlleleDelta 規則
L5   殘餘 PASS FP          11,598       核心問題
     └─ raw_absent         11,430       98.6%，需 Caller level 解決
     └─ persistent         87           0.75%，hardest FP
     └─ raw_filter         63           0.54%，paired filter
     └─ longphase_s        18           0.15%，paired phasing
```

**主要結論**：TO FP 問題的本質是 **ClairS-TO calling model 在無 normal 樣本下的系統性假陽性**，不是甲基化可解決的信號問題。ISM 在此 track 的最大貢獻仍是 FN rescue（H012: +0.009365），而非 FP 過濾。

---

*報告由研究迴圈 v1.1.0 + 手動 FP provenance 分析生成*
*研究者：Claude Sonnet 4.6 + 人類研究員*
*分析腳本：`scripts/analysis/build_to_fp_provenance_analysis.py`*
*輸出數據：`research/fp_provenance/20260322_hcc1395_5khz_to/`*
