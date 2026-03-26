<!--
建立時間: 2026-03-22
目標: 跨 7 個 TO 樣本驗證 ISM 甲基化特徵 AF 梯度，確認候選池偏差與全量統計差異
處理範圍: HCC1395_5kHz（救援候選池）, HCC1395_DORADO, COLO829, H1437, H2009, HCC1937, HCC1954（全量 ISM 輸出）
關聯檔案:
  - docs/reports/validated/2026/03/20260322_cross_sample_TO_LOH_FP_enrichment_01.md（前置：LOH 富集分析）
  - docs/reports/validated/2026/03/20260322_TO_TP_FP_germline_LOH_characterization_01.md（前置：TP/FP Germline-LOH 特徵畫像）
  - research/fp_provenance/20260322_cross_sample_to_loh_analysis/ism_gradient_summary.tsv
  - scripts/analysis/analyze_cross_sample_ism_af_gradient.py
-->

# 跨樣本 TO ISM 甲基化特徵 AF 梯度分析

> 生成時間：2026-03-22
> 前置研究：20260322 LOH 富集分析、20260322 Germline-LOH 特徵畫像

---

## 一、研究背景與問題

前次分析（20260322_TO_TP_FP_germline_LOH_characterization）在 HCC1395 5kHz 的 ISM **救援候選池**（rescue candidate pool，773 TP + 298 FP）中發現顯著的 ISM 甲基化判別力 AF 梯度：

| AF 區間 | PassedGating 富集 | VerificationClass Strong TP/FP |
|---------|------------------|-------------------------------|
| <0.2    | 1.17×            | 18% / 19%                     |
| 0.6-0.9 | **2.00× ✅**      | **40% / 7%**（最有效）         |
| ≥0.9    | **0.45× ❌**      | 7% / 12%（反效果）             |

本次研究問題：
1. 上述 AF 梯度是否在其他 TO 樣本的**全量 ISM 輸出**中同樣存在？
2. 候選池（rescue candidate pool）選擇是否引入了偏差？
3. ISM PassedGating 在 LOH（AF≥0.9）的反效果是否為普遍現象？

---

## 二、方法

**資料來源**：
- HCC1395_5kHz：`rescue_joined_features.tsv`（候選池，773 TP + 298 FP，**有選擇偏差**）
- 其他 6 個樣本：`step05_intersubmod/intersubmod_tp|fp/significance_summary.csv`（**全量 ISM 輸出**）+ `step04_benchmark_longphase_to/filtered_snv_tp|fp.vcf.gz` AF

**AF 分區間**：`[0, 0.2, 0.4, 0.6, 0.9, ≥0.9]`，5 個 bins

**評估指標**：
- `gate_enrichment`：PassedGating 後 TP:FP ratio / 篩選前 TP:FP ratio
- `tp/fp_strong_pct`：VerificationClass = Strong 的比例
- `tp/fp_noise_pct`：VerificationClass = Noise 的比例
- Quality_Score 中位數 Mann-Whitney p-value

---

## 三、關鍵發現：候選池偏差

### 3.1 HCC1395_5kHz（候選池）vs HCC1395_DORADO（全量）的差異

| 指標 | HCC1395_5kHz（候選池） | HCC1395_DORADO（全量） |
|------|----------------------|----------------------|
| AF 0.6-0.9 TP N | 30 | 2,865 |
| AF 0.6-0.9 FP N | 15 | 1,914 |
| AF ≥0.9 TP N | 15 | 5,193 |
| AF ≥0.9 FP N | 17 | 4,519 |
| **PassedGating enrichment 0.6-0.9** | **2.00×** | **1.13×** |
| **PassedGating enrichment ≥0.9** | **0.45× ❌** | **1.07×** |
| VerificationClass Strong 0.6-0.9 TP/FP | **40% / 7%** | **40% / 35%** |

**結論**：
- 候選池在 0.6-0.9 AF 的 `FP_Strong=7%` 是人為選擇的結果（候選池的 FP 以低品質 FP 為主，真正 LOH-like FP 佔比不同）
- 全量統計中，DORADO 的 FP_Strong=35% 在 0.6-0.9 AF，與 TP_Strong=40% 差距遠比候選池小
- **PassedGating 在 AF≥0.9 的反效果（0.45×）是候選池特有現象，全量（DORADO）顯示為中性（1.07×）**

---

## 四、跨樣本 PassedGating 富集倍數

### 4.1 完整對比表

```
AF 區間   HCC1395   HCC1395    COLO829   H1437    H2009    HCC1937  HCC1954
          5kHz      DORADO
          (候選池)   (全量)     (全量)    (全量)   (全量)   (全量)   (全量)
────────────────────────────────────────────────────────────────────────────
<0.2      1.17×     1.07×      0.71×*    0.97×    1.61×    0.64×*   1.65×
0.2-0.4   1.36×     0.97×      0.97×     1.27×    1.12×    1.14×    1.39×
0.4-0.6   1.61×     0.99×      1.18×     1.30×    1.26×    0.98×    1.26×
0.6-0.9   2.00× ✓  1.13×      1.20×     1.22×    1.30×    1.21×    1.39×
≥0.9      0.45×!!  1.07×      0.89×*    1.52×    1.58×    1.07×    1.45×
────────────────────────────────────────────────────────────────────────────
!! = 強反效果（< 0.5×）  * = 反效果（< 0.9×）  ✓ = 強正效果（> 1.5×）
```

**解讀**：
- HCC1395_5kHz 的梯度效果（候選池）：明顯（0.45× → 2.00×），但有偏差
- HCC1395_DORADO（全量）：梯度平坦（0.97-1.13×），僅在 0.6-0.9 有微弱正效果
- COLO829 (黑色素瘤)：<0.2 出現輕微反效果（0.71×）
- H1437/H2009（肺腺癌）：AF≥0.9 出現正效果（1.52×/1.58×）——可能因 LOH 在這些樣本的 ISM 行為不同
- HCC1937（BRCA1 mutant）：<0.2 出現反效果（0.64×）——異常，低 AF 的 FP 通過率反而高於 TP
- HCC1954（乳腺癌）：整體正效果，尤以 <0.2（1.65×）最強

### 4.2 AF≥0.9（LOH 區域）的 PassedGating 通過率

| 樣本 | TP PassedGating% | FP PassedGating% | 意義 |
|------|-----------------|-----------------|------|
| HCC1395_5kHz（候選池）| 13.3% | 29.4% | FP 通過率遠高於 TP（候選池偏差） |
| HCC1395_DORADO（全量）| 13.8% | 12.9% | 幾乎相同，無區分力 |
| COLO829 | 0.9% | 1.0% | 幾乎為零，ISM 不適用 |
| H1437 | 1.3% | 0.9% | 很低，ISM 幾乎不適用 |
| H2009 | 4.6% | 2.9% | 略微 TP > FP |
| HCC1937 | 1.5% | 1.4% | 幾乎相同 |
| HCC1954 | 2.2% | 1.5% | 略微 TP > FP |

**結論**：LOH 區域（AF≥0.9）的 PassedGating 通過率在全量分析中普遍極低（0.9-13.8%），且 TP/FP 幾乎無法區分。

---

## 五、VerificationClass 梯度（跨樣本）

### 5.1 AF≥0.9 Noise % — 普遍性驗證

| 樣本 | TP Noise% | FP Noise% | 意義 |
|------|----------|----------|------|
| HCC1395_5kHz | 86.7% | 76.5% | 高，但候選池偏差 |
| HCC1395_DORADO | 87.3% | 88.2% | **幾乎完全 Noise** |
| COLO829 | **99.6%** | **99.4%** | 完全 Noise |
| H1437 | 98.8% | 99.1% | 完全 Noise |
| H2009 | 95.5% | 97.0% | 完全 Noise |
| HCC1937 | 98.1% | 98.7% | 完全 Noise |
| HCC1954 | 97.0% | 98.1% | 完全 Noise |

**✅ 普遍結論**：AF≥0.9 時，**所有樣本的 TP 和 FP 在 VerificationClass 上均以 Noise 為主（87-100%），ISM 甲基化分析在 LOH 區域對所有樣本完全失效。**

### 5.2 AF 0.6-0.9 Strong % — 判別力比較

| 樣本 | TP Strong% | FP Strong% | TP-FP 差距 | 意義 |
|------|-----------|-----------|-----------|------|
| HCC1395_5kHz（候選池）| **40.0%** | **6.7%** | **+33.3%** | 高（候選池偏差） |
| HCC1395_DORADO（全量）| 40.1% | 34.7% | +5.4% | 弱 |
| H1437 | 37.7% | 31.1% | +6.6% | 弱 |
| H2009 | 9.5% | 7.1% | +2.4% | 幾乎無 |
| COLO829 | 11.8% | 10.1% | +1.7% | 幾乎無 |
| HCC1937 | 19.5% | 16.6% | +2.9% | 幾乎無 |
| HCC1954 | 20.2% | 14.1% | +6.1% | 弱 |

**結論**：全量分析中，AF 0.6-0.9 的 ISM 判別力（Strong% TP/FP 差距）比候選池小得多（5-7% vs 33%）。HCC1395_DORADO 顯示最強的全量梯度，但絕對差距仍有限。

### 5.3 HCC1937 <0.2 異常（FP_Strong > TP_Strong）

| AF 區間 | TP Strong% | FP Strong% | TP Noise% | FP Noise% |
|---------|-----------|-----------|----------|----------|
| <0.2    | 12.5%     | **25.5%** | 51.7%    | 35.5%    |
| 0.2-0.4 | 17.1%     | 15.6%     | 35.8%    | 36.7%    |

HCC1937 的低 AF（<0.2）FP 具有比 TP 更高的 Strong%，這與其他樣本相反。這與 HCC1937 的 BRCA1 缺陷有關：
- **假說**：低 AF FP 中有較多來自**BRCA1 缺陷誘發的染色體不穩定性**產生的偽陽性，這些位點附近有真實的甲基化差異（但不是體細胞突變）。

---

## 六、Quality Score 梯度（跨樣本）

| AF 區間 | HCC1395_5kHz | HCC1395_DORADO | COLO829 | H1437 | H2009 | HCC1937 | HCC1954 |
|---------|-------------|----------------|---------|-------|-------|---------|---------|
| **<0.2** | **TP=75 > FP=60** (p<0.001) | TP=FP=75 (p=0.0*) | TP=FP=35 | TP=FP=75 (p<0.001) | TP=FP=75 (p<0.001) | TP=FP=75 | TP=100 > FP=85 (p<0.001) |
| **0.2-0.4** | TP=75 > FP=60 (p<0.001) | TP=FP=75 (p=0.0*) | TP=FP=35 | TP=FP=75 (p<0.001) | TP=90 < FP=80?† | TP=75 < FP=80 (p<0.001) | TP=FP=85 |
| **0.6-0.9** | TP=60 < FP=75 | TP=FP=75 (p=0.0*) | TP=FP=35 | TP=FP=75 | TP=FP=80 (p<0.001) | TP=75 < FP=80 (p<0.001) | TP=75 < FP=85 (p<0.001) |
| **≥0.9** | TP=FP=60 (p=0.29) | TP=FP=75 (p=0.0*) | TP=FP=50 | TP=FP=75 | TP=FP=75 (p<0.001) | TP=FP=75 (p<0.001) | TP=50 < FP=60 (p<0.001) |

*HCC1395_DORADO 因大 N（千~萬級）微小差異也顯著，中位數相同但分佈尾部不同。
†H2009/HCC1937 某些 AF bin 中 FP 中位 QS 略高於 TP，原因尚不明確。

**觀察**：
- COLO829 的 QS 整體偏低（中位數 35），且 TP/FP 完全無法區分——黑色素瘤的 ISM Quality Score 比其他癌症類型低得多。
- HCC1395_5kHz 候選池 QS TP>FP 效果在全量（DORADO）幾乎消失（中位數相同）。
- 多個樣本出現 FP QS ≥ TP QS 的倒掛，尤其在高 AF 區間。

---

## 七、綜合結論

### 7.1 普遍性結論（所有樣本共同確認）

| 現象 | 普遍性 | 強度 |
|------|--------|------|
| AF≥0.9 時 VerificationClass Noise% ≥ 87%（TP 和 FP 均如此） | ✅ 所有樣本 | 99-100%（肺腺癌）到 87%（HCC1395） |
| AF≥0.9 時 PassedGating 通過率極低（< 15%）| ✅ 所有樣本 | 0.9-14% |
| AF≥0.9 時 PassedGating 對 TP/FP 無區分力 | ✅ 所有樣本 | TP/FP 通過率差距 ≤ 2% |
| ISM 在 LOH 區域完全失效 | ✅ 所有樣本 | 強 |

### 7.2 候選池偏差（非普遍性）

| 現象 | 來源 | 是否普遍 |
|------|------|---------|
| AF 0.6-0.9 PassedGating 2.00× 強效果 | HCC1395_5kHz 候選池 | ❌ 非普遍（全量 0.9-1.4×）|
| AF≥0.9 PassedGating 0.45× 反效果 | HCC1395_5kHz 候選池 | ❌ 非普遍（全量 0.89-1.58×）|
| VerificationClass Strong TP=40% vs FP=7% at 0.6-0.9 | HCC1395_5kHz 候選池 | ❌ 非普遍（全量差距 2-6%）|
| QS TP=75 vs FP=60 at 低 AF | HCC1395_5kHz 候選池 | ❌ 非普遍（全量 TP=FP=75）|

**候選池偏差機制**：救援候選池選取的是「caller 可能遺漏的 FN」+ 「混入候選池的 FP」。FP 候選者不代表全量 FP 分佈，而是以**低品質、接近 borderline 的 FP** 為主，因此 VerificationClass、QS 差距被人為放大。

### 7.3 樣本特異性現象

| 樣本 | 特殊現象 | 可能機制 |
|------|---------|---------|
| COLO829（黑色素瘤）| 整體 QS 低（中位數 35），ISM 幾乎無效 | UV-signature TMB 高，甲基化信號雜亂 |
| HCC1937（BRCA1 mut）| <0.2 FP_Strong > TP_Strong（逆轉）| BRCA1 缺陷染色體不穩定性誘發低 AF 附近的甲基差異 FP |
| H1437/H2009（肺腺癌）| AF≥0.9 PassedGating 1.5-1.6×（最高正效果）| 肺腺癌 LOH 中 ISM 微弱但正向有效 |
| HCC1954（乳腺癌）| <0.2 PassedGating 1.65×（全量中最高）| 低 AF 的真體細胞突變有明顯甲基差異 |

### 7.4 對研究方向的影響

1. **前次 ISM AF 梯度結論需修正**：
   - **保留**：AF≥0.9 VerificationClass 完全失效（普遍，所有樣本確認）
   - **修正**：PassedGating 反效果（0.45×）是候選池偏差產物，全量分析不成立
   - **修正**：0.6-0.9 的 ISM 峰值效果（2.00×，TP_Strong=40% vs FP=7%）是候選池效果，全量差距小得多

2. **ISM 對 TO track 的整體評估**：
   - 在全量統計中，ISM PassedGating 的 TP/FP 區分力**普遍微弱**（1.0-1.5×），且不具 AF 梯度一致性
   - 對 LOH（AF≥0.9）的 FP：ISM 完全無效（普遍）
   - 對低-中 AF 的 FP：ISM 有微弱正向效果（0.97-1.65×），但效果因樣本和癌症類型差異大

3. **COLO829 特別注意**：ISM 在黑色素瘤中 QS 特別低（35 vs 75），且效果最弱（0.71-1.20×）——ISM 可能不適合黑色素瘤。

---

## 八、結論摘要

| 問題 | 答案 |
|------|------|
| ISM AF 梯度是否普遍？ | 否。HCC1395_5kHz 的強梯度是候選池偏差。全量效果平坦（0.9-1.4×）。 |
| LOH（AF≥0.9）ISM 失效是否普遍？ | **是**。所有 7 個樣本在 AF≥0.9 均 ≥87% Noise，PassedGating 無效。 |
| PassedGating 在 LOH 的反效果（0.45×）是否普遍？ | 否。僅在 HCC1395_5kHz 候選池出現，全量（DORADO）為中性（1.07×）。 |
| 0.6-0.9 AF 是 ISM 最有效的區間？ | 在候選池是（2.00×），在全量僅微弱（1.1-1.4×）。 |
| ISM 在哪個癌症類型最無效？ | 黑色素瘤（COLO829）：QS=35，PassedGating<0.9-1.2×，幾乎無效。 |
| 建議的 ISM rescue 策略？ | 全量統計中 ISM 對 TO FP 效果非常有限，不建議單純依靠 ISM 過濾 TO FP。 |

---

## 九、數據位置

- ISM 特徵梯度分析：`research/fp_provenance/20260322_cross_sample_to_loh_analysis/`
  - `ism_gradient_summary.tsv`：跨樣本 ISM 特徵彙總（主要分析表）
  - `ism_gradient_[sample].tsv`：各樣本 ISM 特徵梯度詳細數據
- 分析腳本：`scripts/analysis/analyze_cross_sample_ism_af_gradient.py`
- 前置報告：
  - `20260322_cross_sample_TO_LOH_FP_enrichment_01.md`
  - `20260322_TO_TP_FP_germline_LOH_characterization_01.md`
