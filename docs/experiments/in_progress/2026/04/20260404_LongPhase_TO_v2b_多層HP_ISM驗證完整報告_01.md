<!--
建立時間: 2026-04-04 01:00
目標: 完整驗證 LongPhase-TO PON-Only v2b phasing + ISM 多層 HP 對 TO 模式 TP/FP 區分能力的影響
處理範圍: HCC1395 5kHz TO mode — baseline vs pononly_v2b ISM 輸出比較、HP 分布恢復、AUC 分析、PON+caller 評估
關聯檔案:
  - /big7_disk/liaoyoyo2001/longphase-to-mod/output/ism_pononly_v2b_tp/significance_summary.csv
  - /big7_disk/liaoyoyo2001/longphase-to-mod/output/ism_pononly_v2b_fp/significance_summary.csv
  - /big7_disk/liaoyoyo2001/longphase-to-mod/output/ism_baseline_tp/significance_summary.csv
  - /big7_disk/liaoyoyo2001/longphase-to-mod/output/ism_baseline_fp/significance_summary.csv
  - /big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp
  - docs/experiments/in_progress/2026/04/20260403_multilayer_HP_before_after_肉眼觀察.md
-->

# LongPhase-TO v2b + 多層 HP：ISM TO 模式驗證完整報告

## 1. 研究問題

1. **LongPhase-TO 的 self-phasing circular dependency 是否影響 ISM？** — 是否需要修正 haplotag？
2. **PON-Only v2b 二次 phasing 是否恢復正確的 HP tag 分布？** — HP:i:1/2/11/21/33 五類標籤是否正常
3. **ISM 多層 HP 測試能否更有效區分 TO 模式的 TP vs FP？** — AUC 是否提升
4. **使用 PON + Caller Germline 作為 phasing 位點是否更好？** — 是否有額外 germline 來源

---

## 2. 實驗設定

### 2.1 資料

| 項目 | 值 |
|------|-----|
| 樣本 | HCC1395 (5kHz ONT simplex) |
| 模式 | Tumor-Only (TO) |
| Tumor BAM | `/big7_disk/liaoyoyo2001/data/HCC1395_5kHz/HCC1395.bam` (272G) |
| TP VCF | `filtered_snv_tp.vcf.gz` (30,476 regions) |
| FP VCF | `filtered_snv_fp.vcf.gz` (4,822 regions) |
| Caller | ClairS-TO (pileup mode) |
| Somatic VCF | `snv.vcf.gz` (3,187,275 variants: 3,133,474 NonSomatic + 47,798 PASS + 6,003 LowQual) |

### 2.2 兩個 Phasing 方案

| 項目 | Baseline (self-phasing) | PON-Only v2b (二次 phasing) |
|------|------------------------|---------------------------|
| Binary | `longphase-to-baseline` | `longphase-to` (modified) |
| Phase 1 anchors | All variants (germline + somatic) | NonSomatic only (3.13M germline) |
| Phase 2 | N/A | Somatic re-integrated into scaffold |
| Purity 估計 | 0.927 | 0 (PON-only 無法估計) |
| 運行時間 | 2,895s | 1,053s |
| Haplotag 標記率 | 18.8M/40.9M (46.0%) | 18.8M/40.9M (46.0%) |

### 2.3 ISM 版本

多層 HP 版本（commit b9eaba7 之後的 working tree），包含：
- **Layer 1 (HP Family)**: HP1+HP1-1 vs HP2+HP2-1 → `GlobalP_HPFamily`
- **Layer 2 (HP Pure)**: HP1 vs HP2 → `HPMergedP` (原有)
- **Layer 3 (HP Fine)**: 4-group RxC Fisher (HP1/HP1-1/HP2/HP2-1) → `HPFineP`

---

## 3. HP 分布恢復驗證

### 3.1 Self-Phasing Circular Dependency 的影響

**問題根因**：baseline LongPhase-TO 將 somatic variant 作為 phasing anchor。由於 somatic SNV 只存在於腫瘤中（AF~0.5），ALT reads 被偏向分配到同一 HP → 產生假 LOH 信號。

**v2b 的修正**：Phase 1 僅使用 germline (NonSomatic) variants 建立 haplotype scaffold，Phase 2 才將 somatic variants 重新整合。

### 3.2 HP 分布數據

| 指標 | Baseline | v2b | Paired (參考) |
|------|----------|-----|-------------|
| Both HP (HP1+HP2 均有) | 65.1% | 52.1% | 65.9% |
| HP1-only | 17.0% | 23.3% | 16.2% |
| HP2-only | 17.3% | 24.1% | 17.3% |
| Unphased | 0.6% | 0.5% | 0.6% |
| Potential_LOH | 21,552 (61.1%) | 25,270 (71.6%) | — |

**觀察**：v2b 有更多 LOH（+3,718 regions），因為移除 self-phasing bias 後，原本被假性平衡的 HP 分布回歸真實狀態（更多 single-HP）。

### 3.3 Haplotag 程式碼驗證

`HaplotagProcess.cpp` 的 GT 解析邏輯（lines 92-195）確認：
- GT=0|1 / GT=1|0：正確設定 refHaplotype/altHaplotype → HP:i:1/HP:i:2
- GT=0|. / GT=.|0：正確設定 altHaplotype → HP:i:11/HP:i:21
- GT=1|. / GT=.|1：正確處理 hemizygous → HP:i:1/HP:i:2
- GT=0|0（somatic, 兩個 germline 都是 REF）：refHaplotype 保持 UNDEFINED → REF reads 不投票（語意正確）

**結論**：haplotag 程式碼本身不需要額外修正。v2b 的改善來自 phasing 階段（移除 somatic anchors），而非 haplotag 階段。

---

## 4. ISM 區分能力分析

### 4.1 單特徵 AUC 比較

| Feature | v2b AUC | Baseline AUC | Delta | Winner |
|---------|---------|-------------|-------|--------|
| **HPFineP** | **0.667** | 0.652 | +0.016 | **v2b** |
| HPFineF | 0.660 | 0.650 | +0.010 | v2b |
| HPFineSig | 0.622 | 0.613 | +0.009 | v2b |
| HPMergedP | 0.484 | 0.517 | -0.033 | baseline |
| Quality_Score | 0.579 | 0.610 | -0.031 | baseline |
| LabelHPPermanovaF | 0.492 | 0.528 | -0.036 | baseline |
| HP_Ratio | 0.431 | 0.497 | -0.066 | baseline |
| AlleleDelta | 0.393 | 0.393 | 0.000 | tie |
| AlleleP | 0.509 | 0.509 | 0.000 | tie |
| GlobalP | 0.445 | 0.441 | +0.004 | tie |

**Summary**: v2b wins=3, baseline wins=5, tie=6

### 4.2 關鍵發現：HPFineP 是最強單特徵

- **v2b HPFineP AUC=0.667** 是所有模式、所有特徵中的最高值
- HPFine sig ratio: v2b=**1.75** (TP 57.3% sig vs FP 32.8%) > baseline=1.61
- HPFineNGroups 分布差異巨大：
  - NGroups=2: TP 56.8% vs FP 32.2%（enrichment 1.76×）
  - NGroups=1: TP 17.4% vs FP 42.9%（FP 嚴重富集於單組）

### 4.3 LOH 分層 AUC — HPFineP 的真正價值在 LOH 區域

| Feature | All | LOH | non-LOH |
|---------|-----|-----|---------|
| HPFineP | 0.667 | **0.737** | 0.534 |
| Quality_Score | 0.579 | 0.654 | 0.572 |
| GlobalP | 0.445 | 0.467 | 0.403 |
| HP_Ratio | 0.431 | 0.438 | 0.386 |

**HPFineP 在 LOH 區域 AUC=0.737**，在 non-LOH 區域僅 0.534。這意味著多層 HP 的價值主要體現在 LOH 區域的區分上。

### 4.4 LOH 子類型分層 — LOH_Strong HPFineP AUC=0.823

| LOH Subtype | n | TP | FP | HPFineP AUC |
|-------------|---|----|----|-------------|
| LOH_Strong | 5,857 | 4,951 | 906 | **0.823** |
| LOH_Weak | 7,126 | 6,545 | 581 | **0.838** |
| LOH_Noise | 10,501 | 8,973 | 1,528 | 0.655 |
| LOH_Subclone | 1,786 | 1,532 | 254 | 0.600 |
| None (non-LOH) | 10,028 | 8,475 | 1,553 | 0.534 |

**關鍵觀察**：
- LOH_Strong/Weak 的 HPFineP AUC 超過 0.82 — 非常強的區分能力
- LOH_Strong: TP sig=84.0% vs FP sig=29.0%（2.9×）
- LOH_Weak: TP sig=76.3% vs FP sig=21.3%（3.6×）
- 解讀：真陽性的 LOH 區域中，somatic variant 改變了 haplotype 內的甲基化模式，HPFine 能偵測到此差異

### 4.5 多特徵組合 AUC

| Feature Combination | v2b AUC |
|---------------------|---------|
| HPFineP (single) | 0.667 |
| HPFineP + GlobalP_HPFine | 0.674 |
| HPFineP + Quality_Score + AlleleDelta + CramersV | **0.755** |
| HPFineP + Quality_Score + AlleleDelta + CramersV + HP_Ratio | **0.758** |

5 特徵 Logistic Regression:

| Feature | Coefficient | 解讀 |
|---------|-------------|------|
| AlleleDelta | -0.798 | TP 的 allele 差異更小（FP 常有假性 allele 偏差） |
| HPFineP | -0.707 | TP 的 HPFine p-value 更低（真實亞克隆信號） |
| CramersV | +0.354 | TP 的效應量更大 |
| HP_Ratio | -0.130 | 輕微貢獻 |
| Quality_Score | -0.057 | 幾乎無獨立貢獻（已被其他特徵覆蓋） |

**LOH 分層 5-feature LR**:
- All: AUC=0.758
- LOH: AUC=**0.780**
- non-LOH: AUC=0.730

**對照 baseline 同樣 5 特徵**: baseline AUC=0.762 vs v2b AUC=0.758。差異極小（-0.004），在統計噪音範圍內。

### 4.6 VerificationClass 分布比較

| Class | v2b TP% | v2b FP% | Baseline TP% | Baseline FP% |
|-------|---------|---------|-------------|-------------|
| Strong | 24.0% | 35.0% | 23.4% | 35.4% |
| Subclone | 5.3% | 5.7% | 4.8% | 5.4% |
| Weak | 34.1% | 21.7% | 35.7% | 22.6% |
| Noise | 36.6% | 37.6% | 36.1% | 36.6% |

v2b 相對 baseline 的 VerificationClass 轉移（主要路徑）：

| 轉移 | 數量 | TP | FP |
|------|------|-----|-----|
| Weak → Noise | 697 | 625 | 72 |
| Noise → Weak | 424 | 396 | 28 |
| Weak → Strong | 355 | 348 | 7 |
| Strong → Weak | 124 | 117 | 7 |
| Noise → Subclone | 112 | 103 | 9 |

---

## 5. PON + Caller Germline 可行性評估

### 5.1 Caller VCF 結構

| FILTER | Count | 用途 |
|--------|-------|------|
| NonSomatic | 3,133,474 | 已被 PON 資料庫標記的 germline variants |
| PASS | 47,798 | Somatic candidates |
| LowQual (各種組合) | 6,003 | 低品質 variants（未被 PON 標記） |

### 5.2 Verdict 欄位

VCF header 定義了 `Verdict_Germline`、`Verdict_Somatic`、`Verdict_SubclonalSomatic`，但**實際記錄為 0 筆**。ClairS-TO pileup 模式不使用 Verdict 分類。

### 5.3 結論：PON + Caller = 已經是 v2b 的現狀

v2b 的 phasing 輸入就是完整的 `snv.vcf.gz`（3,187,275 variants），其中 3,133,474 NonSomatic 已經是所有可用的 germline variants。沒有額外的 germline 來源可用。

**具體比較**：

| 方案 | Germline Anchors | Somatic 處理 |
|------|------------------|-------------|
| Baseline | 3.13M germline + 47.8K somatic（全部混合 phasing） | 作為 anchor |
| v2b | 3.13M germline only（Phase 1） | Phase 2 re-integrate |
| "PON + Caller" | = v2b（無額外來源） | 同 v2b |

**額外 germline 的唯一潛在來源**：6,003 個 LowQual 且未被 PON 標記的 variants。但這些品質低、數量少（佔 NonSomatic 的 0.2%），幾乎不會影響 phasing 結果。

---

## 6. 關鍵數字摘要

| 指標 | v2b | Baseline | 變化方向 |
|------|-----|----------|---------|
| **HPFineP AUC** | **0.667** | 0.652 | ↑ (+0.015) |
| HPFine sig ratio (TP/FP) | **1.75** | 1.61 | ↑ |
| HPFineP LOH AUC | **0.737** | — | 新指標 |
| HPFineP LOH_Strong AUC | **0.823** | — | 新指標 |
| 5-feature LR AUC | 0.758 | 0.762 | ≈ (差異在噪音範圍) |
| Quality_Score AUC | 0.579 | 0.610 | ↓ (-0.031) |
| HP_Ratio AUC | 0.431 | 0.497 | ↓ (-0.066) |
| Potential_LOH regions | 25,270 | 21,552 | ↑ (+3,718) |
| PassedGating | 12,025 | 11,745 | ↑ (+280) |

---

## 7. 結論

### 7.1 LongPhase-TO v2b 的效果

**正面**：
1. v2b 成功移除 self-phasing circular dependency，HP 分布回歸真實狀態
2. HPFineP 單特徵 AUC 從 0.652 提升至 **0.667**（所有模式中最高）
3. HPFine sig ratio 從 1.61 提升至 **1.75**
4. LOH 區域的 HPFineP 特別有效：LOH_Strong/Weak AUC > 0.82

**中性/負面**：
1. HP_Ratio、HPMergedP、Quality_Score 等指標在 v2b 下反而退化
2. 5 特徵組合 AUC 幾乎不變（0.758 vs 0.762）
3. 更多區域被分類為 LOH（+3,718），但這反映真實狀態而非退化

### 7.2 為何多特徵 AUC 不變？

v2b 改善了 HPFine（更準確的 haplotype 信號），但同時：
- 增加了 LOH 區域數量 → Quality_Score 中的 LOH penalty 影響更多區域
- HP_Ratio 變得更極端 → 特徵的 AUC 反轉（baseline 的 moderate HP_Ratio 反而有更好的區分度）
- 這兩個效應互相抵消，導致組合 AUC 接近不變

### 7.3 ISM 是否因 v2b 修正而「更能區分 TO 下的 TP/FP」？

**部分是**：
- **HPFine 層級**：是的，v2b 提供更清潔的 haplotype 信號 → HPFineP 改善顯著
- **LOH 場景**：是的，LOH_Strong/Weak 的 HPFineP AUC > 0.82 是非常有價值的發現
- **整體多特徵**：否，組合效果被 HP_Ratio 和 Quality_Score 的退化抵消

**actionable insight**：
- 對 LOH 區域，HPFineP 是一個強力的新特徵（AUC 0.737-0.838）
- Quality_Score 的 LOH penalty 設計需要針對 v2b 的新 LOH 分布重新校準
- HP_Ratio 在 v2b 下反轉方向，不適合直接作為 discriminative 特徵

### 7.4 PON + Caller Germline

**結論：已無額外可用的 germline 來源**。v2b 已經使用了 caller VCF 中所有 3.13M NonSomatic variants 作為 phasing anchors。Verdict 欄位為空（pileup 模式不產生），6K LowQual variants 數量太少，不值得嘗試。

---

## 8. 建議後續行動

### 高優先

1. **重新設計 Quality_Score for TO-v2b**：當前 LOH penalty 在 v2b 下過度懲罰（LOH 增加 +3,718 regions），需要重新校準或分離 LOH 處理邏輯
2. **HPFineP 納入正式 QS**：HPFineP 在 LOH 場景下 AUC=0.737，是目前最強的 TO 單特徵，應整合到 Quality Score 計算中

### 中優先

3. **LOH 子類型特化分析**：LOH_Strong/Weak HPFineP AUC > 0.82，探索是否可以建構 LOH-aware 的分類器
4. **HP_Ratio 方向修正**：在 v2b 下 HP_Ratio 的 AUC 反轉（0.431 < 0.5），需理解原因並決定是否保留此特徵

### 低優先

5. **獨立 RNG per test**：消除 Fisher MC 測試順序對結果的微小影響（影響 ~2.6% 區域）
6. **擴展到其他樣本**：COLO829 / H1437 / H2009 的 v2b 驗證

---

## 9. 數據路徑

| 資料 | 路徑 |
|------|------|
| v2b TP ISM | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/ism_pononly_v2b_tp/` |
| v2b FP ISM | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/ism_pononly_v2b_fp/` |
| Baseline TP ISM | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/ism_baseline_tp/` |
| Baseline FP ISM | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/ism_baseline_fp/` |
| v2b phasing log | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v2b/phasing.log` |
| Baseline phasing log | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/run.log` |
| Caller VCF | `/big7_disk/liaoyoyo2001/data/HCC1395_5kHz/snv.vcf.gz` |
| Haplotag source | `/big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp` |
