<!--
建立時間: 2026-03-22
目標: TO 全量 TP 特徵畫像 + FP Germline/LOH/Artifact 分類 + ISM 甲基化特徵差異
處理範圍: HCC1395 5kHz TO - TP(28,509) + FP(11,606), ISM 候選池(773 TP + 298 FP)
關聯檔案:
  - docs/reports/validated/2026/03/20260322_TO_FP_provenance_methylation_analysis_01.md（前置研究）
  - research/fp_provenance/20260322_hcc1395_5khz_to/
  - step04_benchmark_longphase_to/tp.vcf / fp.vcf
  - rescue_joined_features.tsv（ISM 候選池特徵）
-->

# TO TP 特徵畫像與 FP Germline/LOH/Artifact 分類研究

> 生成時間：2026-03-22
> 資料集：HCC1395 5kHz TO (5kHz_simplex_5mCG_5hmCG, 100% purity cell line)
> 前置研究：FP Provenance Analysis（raw_absent=98.6%，GQ/QS 過濾不可行）

---

## 一、研究背景

前次 FP Provenance 分析（同日）確認：
- TO 殘餘 PASS FP = 11,598，其中 98.6%（11,430 個）為 raw_absent（ClairS-TO 特有 calls）
- GQ 和 Quality_Score 分佈在 TP/FP 之間完全相同，**任何過濾規則均損失 TP**
- TO FP 過濾方向已耗盡

本研究延伸探索：
1. **TO TP 全量特徵**：AF 分佈、Subclonal/Somatic/LOH 亞型分類
2. **TO FP 的 Germline/LOH/Artifact 性質**：為何 98.6% 的 FP 是 paired mode 從未產生的 calls？
3. **ISM 甲基化特徵**：不同亞型的 FP 是否有可區分的甲基化訊號？
4. **方法論**：能否用 ISM 特徵識別 LOH-germline FP？

**HCC1395 關鍵背景**：
- 100% 純腫瘤細胞系（TNBC, ductal carcinoma）
- 已知有廣泛 LOH（大面積 TP53 mutation + CN alteration）
- 配對 normal = HCC1395BL（同患者 B-lymphoblast）

---

## 二、TO TP 全量特徵畫像

### 2.1 AF 分佈與亞型分類（N=28,509）

| 亞型 | AF 範圍 | Count | % | 生物學解讀 |
|------|---------|-------|---|-----------|
| Subclonal TP | < 0.2 | 2,680 | **9.4%** | 亞克隆體細胞突變（低頻率） |
| Somatic TP | 0.2–0.6 | 17,931 | **62.9%** | 主克隆體細胞突變（主體） |
| Mid-High TP | 0.6–0.9 | 2,917 | **10.2%** | CN gain 或高克隆性突變 |
| LOH-like TP | ≥ 0.9 | 4,981 | **17.5%** | LOH 區域體細胞突變（alt allele 固定） |

**AF 統計**：min=0.07, Q1=0.286, median=0.427, Q3=0.634, max=1.0

```
TP AF 分佈（示意）
   18% ┤                                  ████
   16% ┤     ████
   14% ┤     ████  ████
   12% ┤     ████  ████  ████
   10% ┤     ████  ████  ████  ████
    8% ┤     ████  ████  ████  ████  ████
    0%  ──────────────────────────────────
         0.1  0.2  0.3  0.4  0.5  ≥0.9
```

**解讀**：HCC1395 作為 100% 純腫瘤 cell line，62.9% TP 落在 0.2-0.6 range（典型 clonal heterozygous somatic mutation）。**17.5% TP 有 AF ≥ 0.9**，這些是 LOH region 中的體細胞突變——LOH 導致 normal allele 缺失，僅剩含突變的 alt allele，因此 AF 趨近 1.0。這是 HCC1395 廣泛 LOH 的直接體現。

### 2.2 GQ 分佈

| GQ 閾值 | TP 通過 | TP% |
|---------|---------|-----|
| GQ ≥ 3 | 28,509 | 100.0% |
| GQ ≥ 10 | 27,653 | 97.0% |
| GQ ≥ 20 | 12,712 | 44.6% |

**LOH-like TP 的 GQ 特性**：GQ ≥ 20 = 64.8%（高於整體 TP 的 44.6%）。LOH 區域的突變因 alt allele 完全取代，ClairS-TO 讀取到高度一致的 alt reads，因此 GQ 普遍偏高。

### 2.3 H Flag（單一 haplotype 支持）

TP 中 H flag = 58.6%（16,695/28,509）。H flag 意義：ClairS-TO 判斷該 variant 主要支撐在單一 HP（haplotype）上。在 TO 下，haplotagging 不均衡，H flag 在 TP/FP 之間無判別力（FP H flag = 61.9%，與 TP 無顯著差異）。

### 2.4 ISM 候選子集 TP（N=773）

來自 rescue_joined_features.tsv 的 773 個 TP（ISM 已分析的 FN 候選）：

| 亞型 | N | QS 中位數 | PMD 中位數 | CV 中位數 | Noise% |
|------|---|----------|------------|----------|--------|
| **TP (全量)** | **773** | **75.0** | **0.210** | **0.000** | **34.7%** |
| TP Somatic(0.2-0.6) | 393 | 75.0 | 0.208 | 0.000 | 29.5% |
| TP LOH-like(≥0.9) | 15 | 60.0 | 0.178 | 0.000 | 86.7% |

注意：ISM 候選 TP 的 AF < 0.6 佔多數（773 中 Subclonal+Somatic = 728, 94.2%），這是因為 ISM rescue 候選主要是低 GQ 的 FN（被 H012 rescue 的目標）。

---

## 三、TO 殘餘 FP 的 Germline/LOH/Artifact 分類

### 3.1 核心發現：AF 分佈的雙峰結構

| 亞型 | AF 範圍 | FP Count | FP% | TP 對照 |
|------|---------|----------|-----|---------|
| Subclonal-like | < 0.2 | 444 | **3.8%** | 9.4% |
| Somatic-like | 0.2–0.6 | 4,797 | **41.3%** | 62.9% |
| Mid-High | 0.6–0.9 | 2,011 | **17.3%** | 10.2% |
| **LOH-like** | **≥ 0.9** | **4,354** | **37.5%** ← | 17.5% |

**FP AF 統計**：min=0.07, Q1=0.438, median=0.650, Q3=0.967, max=1.0

**關鍵觀察**：
- FP 中 LOH-like（AF≥0.9）佔 **37.5%**，遠高於 TP 的 17.5%（**2.14 倍**）
- FP 的 AF 中位數 = 0.650，明顯高於 TP 的 0.427
- FP GQ 中位數 = 19.0，與 TP 的 19.0 完全相同

### 3.2 LOH Germline-Escape 假說（核心解釋）

**假說**：殘餘 FP 中高比例（~37.5%）是「LOH 區域中的 germline 雜合 SNP」：

```
生物學機制：
正常細胞：germline het SNP → ref/alt 各一 → AF ≈ 0.5
   ↓ 腫瘤演化（HCC1395 廣泛 LOH）
腫瘤細胞：LOH 發生 → ref allele 缺失 → 僅剩 alt allele → AF ≈ 1.0

ClairS-TO 看到：AF=1.0（高信心）→ 誤判為體細胞突變
ClairS-paired 看到：normal AF=0.5 → 立即識別為 germline → 不輸出
Truth set (SEQC2)：這不是體細胞突變 → 標記為 FP
```

**為什麼 paired_raw_naf = 0（raw_absent）**：
- paired mode ClairS 在看到 normal AF=0.5 後，根本不生成這個位點的體細胞 VCF 記錄
- 所以 paired VCF 中完全沒有這個 variant → `paired_raw_naf = 0`（無 record，非 no alt reads）

**支持證據**：
1. FP LOH-like 比例（37.5%）是 TP LOH-like（17.5%）的 2.14 倍
2. 所有 raw_absent FP 的 ndp=0（paired 完全沒有記錄）
3. FP GQ 在 LOH-like 亞型中也很高（GQ≥20 = 60.6%）——因為 alt reads 高度一致
4. HCC1395 是已知廣泛 LOH 的 TNBC 細胞系

### 3.3 FP 亞型特徵比較

| FP 類型 | Count | AF% | GQ 中位數 | GQ≥20% | 描述 |
|---------|-------|-----|----------|--------|------|
| LOH-like (AF≥0.9) | 4,354 | 37.5% | 20.0 | **60.6%** | Germline-in-LOH（高信心偽陽） |
| Somatic-like (0.2-0.6) | 4,797 | 41.3% | 19.0 | 39.6% | 混雜類型（germline het + artifacts） |
| Mid-High (0.6-0.9) | 2,011 | 17.3% | 20.0 | 51.6% | 部分 LOH + CN gain 區域 |
| Subclonal-like (<0.2) | 444 | 3.8% | 13.0 | **11.9%** | 低 AF artifact（最易過濾但數量少）|

**最難過濾的 FP**：LOH-like 和 Mid-High 合計 54.8%，具有高 AF + 高 GQ，ClairS-TO 誤判信心極高。

### 3.4 ClairS-TO Verdict Germline 的局限

ClairS-TO 的 `verdict_germline=True` flag 在 ISM 候選池 FP 中出現率 = **0/298 (0%)**。

這意味著：ClairS-TO 自身**從未識別**這些逃過 PON 的 germline-in-LOH 為 germline。
- 原因：ClairS-TO 缺少 normal sample，無法判斷 AF=1.0 是 somatic clone 還是 LOH germline
- 唯一能識別的方式是對比 normal BAM

---

## 四、ISM 甲基化特徵的判別能力分析

### 4.1 全量 ISM 候選池特徵比較（298 FP vs 773 TP）

| 組別 | N | med_QS | med_PMD | med_CV | Noise% | med_GQ |
|------|---|--------|---------|--------|--------|--------|
| **TP** | 773 | **75.0** | **0.210** | **0.000** | 34.7% | 14.0 |
| **FP (all)** | 298 | 60.0 | 0.172 | 0.000 | 33.2% | 7.0 |
| FP Subclonal(<0.2) | 159 | 60.0 | 0.155 | 0.000 | 30.8% | 6.0 |
| FP Somatic(0.2-0.6) | 107 | 60.0 | 0.173 | 0.000 | 25.2% | 7.0 |
| FP Mid-High(0.6-0.9) | 15 | 75.0 | 0.201 | 0.000 | 66.7% | 18.0 |
| **FP LOH-like(≥0.9)** | **17** | **60.0** | **0.172** | **0.000** | **76.5%** | **20.0** |

### 4.2 Mann-Whitney U 顯著性測試

| 比較 | QS p-value | CV p-value | PMD p-value | GQ p-value |
|------|-----------|-----------|------------|-----------|
| FP_all vs TP | **<0.0001** | 0.0479 | **0.0001** | **<0.0001** |
| **FP_LOH vs TP** | **0.255** | **0.301** | **0.529** | **<0.0001** |
| FP_Somatic vs TP | 0.0009 | 0.573 | 0.0461 | **<0.0001** |
| FP_Subclonal vs TP | <0.0001 | 0.077 | <0.0001 | **<0.0001** |

**關鍵發現**：
- **FP_LOH（AF≥0.9）的 QS、CV、PMD 與 TP 無顯著差異（p>0.25）**
- 唯一顯著差異：GQ（因為 FP_LOH 在 ISM 候選池中 GQ 偏高，p<0.0001，但方向相反於過濾需求）
- FP_Subclonal 和 FP_Somatic 的 QS、PMD 較低（p<0.001），但數量更少的 FP_LOH 才是主體威脅

### 4.3 VerificationClass 分佈差異

| VerificationClass | TP | FP_all | FP_Somatic | FP_LOH-like |
|------------------|----|--------|------------|-------------|
| Noise | 35% | 33% | 25% | **76%** |
| Weak | 26% | 18% | 24% | 0% |
| Subclone | 3% | 2% | 2% | 6% |
| Strong | 24% | 18% | 19% | 12% |

**FP LOH-like 的 VerificationClass = Noise 高達 76%**。

這是因為：
- LOH-germline FP 在 tumor 中 AF=1.0（只有 alt reads）
- ISM 看到的甲基化分析中，兩個 HP 都是 alt allele → 沒有 HP1 vs HP2 甲基差異
- → AlleleDelta ≈ 0, CramersV ≈ 0 → ISM 判斷為 Noise

然而，**Noise 過濾無效**：TP 中也有 34.7% 為 Noise，大規模過濾 Noise 會損失大量 TP。

---

## 五、LOH 區域 TP vs FP 的特徵對比

### 5.1 LOH-like 亞型的 TP vs FP 完整比較

| 指標 | LOH-like TP (AF≥0.9) | LOH-like FP (AF≥0.9) |
|------|---------------------|---------------------|
| Count（全量 VCF） | 4,981（TP 的 17.5%）| 4,354（FP 的 37.5%）|
| GQ 中位數 | 21.0 | 20.0 |
| GQ ≥ 20% | 64.8% | 60.6% |
| ISM 候選子集 N | 15 | 17 |
| ISM QS 中位數 | 60.0 | 60.0 |
| ISM PMD 中位數 | 0.178 | 0.172 |
| ISM VC=Noise% | 86.7% | 76.5% |

**發現**：LOH-like TP 和 LOH-like FP 的特徵**幾乎完全相同**。
- GQ 分佈相同（高信心）
- ISM QS 相同（中等品質）
- ISM VC=Noise% 相似（兩者都高）
- **→ 在 TO 模式下，沒有任何特徵能區分 LOH 區域的體細胞突變（TP）和 LOH 區域的 germline SNP（FP）**

### 5.2 根本性限制：TO 缺乏 Normal 資訊

```
LOH-like TP（真體細胞突變，AF≥0.9）：
  - 正常細胞：ref/ref（純合 ref）
  - 突變：ref → alt（體細胞突變）
  - LOH：ref allele 進一步缺失 → 只剩 alt → AF=1.0
  → 真正的「somatic + LOH」組合

LOH-like FP（germline het SNP，AF≥0.9）：
  - 正常細胞：ref/alt（germline 雜合）
  - LOH：ref allele 缺失 → 只剩 germline alt → AF=1.0
  → 「germline + LOH」組合（不是體細胞突變）

兩者在 TO 下的觀測值完全相同：
  - tumor BAM：alt reads 完全主導，AF=1.0
  - 沒有 normal sample → 無法判斷正常細胞是 ref/ref 還是 ref/alt
  → TO 模式從根本上無法區分這兩類！
```

---

## 六、結論

### 6.1 主要發現總結

| 問題 | 答案 |
|------|------|
| TO FP 的主要來源是什麼？ | **LOH 區域的 germline SNP（37.5%）+ Somatic-like artifact（41.3%）** |
| FP 和 TP 在甲基化特徵上有差異嗎？ | **FP_all 有（QS、PMD 較低），但 FP_LOH 無（p>0.25）** |
| ISM 能識別 Germline-in-LOH FP 嗎？ | **不能**。LOH-like FP VerificationClass = Noise=76%，但 TP 也有 34.7% Noise |
| GQ 能區分 LOH TP 和 LOH FP 嗎？ | **不能**。GQ 分佈幾乎完全相同（均為高信心） |
| 有任何特徵可過濾 LOH FP 嗎？ | **沒有**。唯一解法是 Normal sample（paired mode） |
| ClairS-TO verdict_germline 有幫助嗎？ | **完全無效**（0/298 FP 被標記為 germline）|

### 6.2 各 FP 亞型的過濾可行性

| FP 亞型 | 數量 | 可過濾性 | 說明 |
|---------|------|---------|------|
| **LOH-like (AF≥0.9)** | **4,354 (37.5%)** | ❌ 不可過濾 | 需 normal sample；ISM/GQ 均無法區分 |
| Mid-High (0.6-0.9) | 2,011 (17.3%) | ❌ 不可過濾 | 同上，部分 LOH，部分 CN gain |
| Somatic-like (0.2-0.6) | 4,797 (41.3%) | ⚠️ 極困難 | FP_Somatic QS 略低，但 TP 誤傷率高 |
| Subclonal-like (<0.2) | 444 (3.8%) | ⚠️ 理論可行 | 低 GQ（中位=13），但只佔 3.8%，delta_F1 微小 |

### 6.3 ISM 在 TO Track 的角色重新定位

```
【確認的 ISM 功能】
✅ FN Rescue（最大價值）：
   - GQ>=3 (H012)：rescue 728 FN → delta_F1 = +0.009365
   - ISM 的 methylation 特徵可輔助識別值得 rescue 的 FN

❌ FP 過濾（已確認無效）：
   - LOH-like FP（37.5%）：無甲基特徵差異
   - Somatic-like FP（41.3%）：QS/PMD 稍低但整體 TP 誤傷風險高
   - GQ、H flag、VerificationClass 均無法有效區分
   - ISM 只覆蓋 2.5% FP（298/11,598），覆蓋率本身就是瓶頸

【根本限制】
TO 模式缺乏 Normal sample，而 LOH 驅動的 FP（最大宗 37.5%）需要 Normal sample 才能識別。
這不是 ISM 算法問題，是 TO pipeline 的系統性限制。
```

---

## 七、下一步建議

### 優先 P1：確認 DORADO TO 是否有相同的 LOH FP 模式

HCC1395 DORADO TO 使用相同的 HCC1395 tumor sample（不同測序平台），
若 LOH 分佈相同，則 DORADO TO FP 應有相似的 AF ≥ 0.9 比例。
此確認可強化 LOH Germline-Escape 假說的跨平台一致性。

### 優先 P2：擴大 ISM FN 覆蓋率（最高效率方向）

```
當前：ISM 覆蓋 773/11,051 FN = 7%
目標：覆蓋 3,000-5,500 FN（30-50%）
潛力：若 rescue precision ≈ 74%（H012 水準），額外 TP ≈ 2,220-4,070
預估 delta_F1：+0.03 至 +0.05
```

### 優先 P3：Subclonal-like FP 的低 GQ 過濾探索

FP 中 Subclonal-like（AF<0.2，N=444）的 GQ 中位數=13，GQ≥20 僅 11.9%，
相比 TP Subclonal GQ≥20=19.1%。雖然規模小，但可能有小幅 F1 改善空間，
值得在 ISM 候選池的 Subclonal FP 進行針對性規則測試。

---

## 八、方法論與限制

### 8.1 資料來源與分析方法

| 分析層次 | 資料 | N |
|---------|------|---|
| 全量 VCF 特徵 | tp.vcf + fp.vcf | 28,509 TP + 11,606 FP |
| Normal AF 分析 | provenance_master.tsv.gz | 11,598 FP（ndp=0: 11,430, ndp>0: 168） |
| ISM 甲基化特徵 | rescue_joined_features.tsv | 773 TP + 298 FP |
| 統計方法 | Mann-Whitney U test（雙尾） | — |

### 8.2 限制

1. **ISM 覆蓋率限制**：ISM 只分析 298/11,598 FP（2.5%），且這些是 ISM rescue 研究的特定子集，不代表全量 FP 的甲基化特性
2. **Normal AF 不可用**：98.5%（11,430/11,598）的殘餘 FP 屬 raw_absent，paired 無記錄，無法直接驗證 normal BAM 的 alt 讀取率
3. **LOH 假說是推論**：沒有直接驗證（需查詢 normal BAM 或 germline VCF 確認 germline het SNP 位置），但所有間接證據一致支持
4. **TO haplotagging 限制**：AlleleDelta/CramersV 在 TO 中已知無判別力（H001-H003 前次驗證），本次 ISM 分析中 CV 中位數均=0 已確認
5. **HCC1395 特異性**：HCC1395 是 TNBC cell line，LOH 程度可能高於一般腫瘤樣本，此比例不一定能外推到其他樣本

---

## 附：數值摘要表

### TP vs FP 全量對比

| 指標 | TP (N=28,509) | FP (N=11,606) |
|------|--------------|--------------|
| AF Q1 | 0.286 | 0.438 |
| AF 中位數 | **0.427** | **0.650** |
| AF Q3 | 0.634 | 0.967 |
| GQ 中位數 | 19.0 | 19.0 |
| LOH-like (AF≥0.9) | 17.5% (4,981) | **37.5% (4,354)** |
| Somatic (0.2-0.6) | 62.9% (17,931) | 41.3% (4,797) |
| Subclonal (<0.2) | 9.4% (2,680) | 3.8% (444) |
| H flag | 58.6% | 61.9% |

---

*報告由手動 Python 分析 + 統計測試生成*
*研究者：Claude Sonnet 4.6 + 人類研究員*
*分析腳本：scripts/analysis/analyze_to_tp_fp_characterization.py*
*數據基礎：step04_benchmark_longphase_to/tp.vcf, fp.vcf + rescue_joined_features.tsv*
