<!--
建立時間: 2026-04-04 17:10
目標: 深度驗證 LOH 區域 AlleleDelta 的生物學機制、F1 改善瓶頸分析、ISM LOH vs LOH.bed 比較
處理範圍: v2b TO ISM (HCC1395 5kHz, 30,476 TP + 4,822 FP) + LOH.bed (pononly_v2b)
關聯檔案:
  - output/synthesis/concluded/loh_verification_analysis/ (7 張圖)
  - docs/experiments/in_progress/2026/04/20260404_LOH_Strong_Weak_7feature_AUC驗證報告_01.md
  - src/core/GlobalTest.cpp (test_allele, lines 149-174)
  - src/core/RegionProcessor.cpp (compute_hp_ratio, lines 38-49)
-->

# LOH 區域 AlleleDelta 深度驗證報告

## 1. 關鍵澄清：AlleleDelta 不是 HP tag 差異

### 1.1 AlleleDelta 實際定義

**AlleleDelta = ALT-read 群 vs REF-read 群的甲基化差異**，不是 HP1 vs HP2。

程式碼路徑：
- `GlobalTest::test_allele()` (`GlobalTest.cpp:149-174`)
- 用 `FullLabel::get_allele_label()` 分類：`ALLELE_ALT` vs `ALLELE_REF`
- 每個 read 的 allele 標籤取決於它是否攜帶 somatic SNV 的突變等位基因

| 欄位 | 分類依據 | 比較 |
|------|---------|------|
| **AlleleDelta** | Somatic SNV 的 ALT/REF allele | ALT-reads vs REF-reads 甲基化 |
| HPMergedDelta | HP tag (1 vs 2) | HP1-reads vs HP2-reads 甲基化 |
| HPFineP | HP tag (1, 1-1, 2, 2-1) | 4-group Fisher 測試 |

### 1.2 ISM 不讀取 LOH.bed

ISM 自行計算 LOH：
```
HP_Ratio = HP1FamilyN / (HP1FamilyN + HP2FamilyN)  // + Laplace smoothing
Potential_LOH = HP_Ratio < 0.1 或 HP_Ratio > 0.9
```

LongPhase-TO 的 `tumor_phased_LOH.bed` 是外部參考，ISM 不使用它。

---

## 2. LOH 下 ALT/REF Read Ratio 是否極端偏斜？

### 2.1 HP Read 分布數據

| 指標 | LOH 內 TP | LOH 內 FP |
|------|-----------|-----------|
| Total HP reads (mean) | 60.0 | 41.9 |
| **Minor HP = 0 (完全單 HP)** | **66.0%** | **67.5%** |
| Minor HP ≤ 2 | 90.1% | 91.4% |
| HP_Ratio ≥ 0.95 | 58.1% | 72.0% |
| HP_Ratio ≥ 0.99 | 48.3% | 58.9% |

**LOH 內 ~66% 的區域完全只有一個 HP 的 reads**，~90% 的區域 minor HP ≤ 2 reads。

### 2.2 為何 LOH 內還會有多個 HP tag？

![LOH HP Biology](../../../../../output/synthesis/concluded/loh_verification_analysis/fig6_loh_hp_biology.png)

**生物學解釋：**

LOH (Loss of Heterozygosity) 有幾種類型：

| LOH 類型 | 機制 | HP tag 影響 |
|---------|------|------------|
| **Deletion LOH** | 一條染色體片段刪除 | 幾乎所有 reads → 保留 HP，minor HP ≈ 0 |
| **cn-LOH (UPD)** | 一條拷貝被另一條取代 | 基因組雖然兩拷貝，但序列相同 → phasing 難以區分 → 可能有少量錯標 |
| **亞克隆 LOH** | LOH 只在部分腫瘤細胞中 | 非 LOH 細胞仍保留兩個 HP → minor HP > 0 |
| **正常細胞汙染** | Tumor purity < 100% | 正常細胞保留兩個 HP → minor HP reads 來自正常細胞 |

**HCC1395 purity ~93%** → 約 7% reads 來自正常細胞 → 解釋了 minor HP 0-2 reads 的來源。

### 2.3 那 AlleleDelta 怎麼算？

即使 LOH 內 HP 極度偏斜，**AlleleDelta 不看 HP tag，看 ALT/REF**：
- 每個 read 是否攜帶 somatic SNV 的 ALT allele → 這取決於 variant caller 標記的 SNV position
- 在 LOH 中，如果 somatic variant 存在，VAF 通常很高（ALT >> REF）
- 但少數 REF reads 仍存在（來自正常細胞或亞克隆）

**所以 ALT/REF ratio 在 LOH 中確實極端偏斜，但仍有少數 REF reads 可供比較。**

---

## 3. LOH 下 AlleleDelta 區分 TP/FP 的生物學機制

![AlleleDelta Distribution](../../../../../output/synthesis/concluded/loh_verification_analysis/fig1_allele_delta_distribution.png)

### 3.1 數據確認

| 區域 | TP AlleleDelta (median) | FP AlleleDelta (median) | FP/TP 比值 |
|------|------------------------|------------------------|-----------|
| **LOH 內** | **0.0138** | **0.0264** | **1.9×** |
| LOH 外 | 0.0224 | 0.0732 | 3.3× |

**LOH 內和 LOH 外，FP 的 AlleleDelta 都比 TP 大。** 差異在 LOH 外更明顯。

### 3.2 生物學解釋

| 類型 | LOH 下的狀況 | ALT/REF 與 Haplotype 的關係 | AlleleDelta |
|------|-------------|---------------------------|-------------|
| **FP germline** | 雜合 germline SNP 在 LOH 中 → 一個 allele 的拷貝數增加 | ALT = 一條 haplotype，REF = 另一條 haplotype → **代表兩條不同的表觀遺傳編程** | **大**（真正的 ASM） |
| **TP somatic** | 體細胞突變在 LOH 後的單一 allele 上 | ALT = 帶突變的 read，REF = 不帶突變的 read → **同一條 haplotype 的兩種 read** | **小**（無天然 ASM） |

**核心差異：**
- **FP germline 的 ALT vs REF 代表兩條真正的 haplotype** → 天然存在 allele-specific methylation (ASM)
- **TP somatic 的 ALT vs REF 不代表 haplotype 差異** → 只是同一背景下的隨機波動

### 3.3 這是「甲基化新訊號」嗎？

**部分是，部分不是：**

| 面向 | 判斷 | 說明 |
|------|------|------|
| 甲基化差異確實存在 | ✅ 是 | FP germline 的 ASM 是真實的生物學現象 |
| 是 ISM 獨有的發現 | ❌ 否 | ASM 與 VAF/AF 高度相關，AF 資訊本身就能提供類似區分 |
| 超越 AF 的額外資訊 | ⚠️ 微弱 | O12 研究顯示控制 AF 後殘餘信號消失 |
| 可推廣到其他場景 | ⚠️ 未知 | 低純度樣本中 minor HP reads 更多，但 ASM 信號也更弱 |

---

## 4. F1 改善為何很小 — 絕對數字分解

![F1 Breakdown](../../../../../output/synthesis/concluded/loh_verification_analysis/fig4_f1_breakdown.png)

### 4.1 基礎數字

| 項目 | 數值 |
|------|------|
| Total TP | 30,476 |
| Total FP | 4,822 |
| FP 佔比 | 13.7% |
| Baseline F1 | **0.9267** |
| 理論完美 F1 | 1.0000 |
| **最大改善空間** | **0.0733** |

### 4.2 LOH_Strong+Weak 覆蓋

| 項目 | 數值 |
|------|------|
| LOH_S+W 內 TP | 11,496 |
| LOH_S+W 內 FP | 1,487 |
| **FP 覆蓋率** | **30.8%** (1,487 / 4,822) |

### 4.3 不同場景的 F1

| 場景 | FP 移除 | TP 損失 | Precision | Recall | F1 | ΔF1 |
|------|---------|---------|-----------|--------|-----|-----|
| Baseline | 0 | 0 | 0.8634 | 1.0000 | 0.9267 | — |
| **LOH_S+W optimal** | **929** | **390** | **0.8854** | **0.9872** | **0.9336** | **+0.007** |
| LOH_S+W aggressive | 1,189 | 1,149 | 0.8898 | 0.9623 | 0.9246 | **-0.002** |
| 完美移除所有 LOH_S+W FP | 1,487 | 1,724 | 0.8961 | 0.9434 | 0.9191 | -0.008 |
| **理論：移除所有 FP** | **4,822** | **0** | **1.0000** | **1.0000** | **1.0000** | **+0.073** |

### 4.4 瓶頸分析

**F1 改善很小的兩個根因：**

| 根因 | 貢獻 | 說明 |
|------|------|------|
| **① FP 覆蓋率低** | 主因 | LOH_S+W 只覆蓋 30.8% FP (1,487/4,822)。即使完美過濾，最多只改善 0.0214 |
| **② TP 連帶損失** | 次因 | optimal threshold 移除 929 FP 的同時損失 390 TP → 淨效果被 TP 損失抵消 |

**效率分析：**
- LOH_S+W 完美過濾上限：F1 = 0.9481 (ΔF1 = +0.0214)
- 實際 optimal 結果：F1 = 0.9336 (ΔF1 = +0.007)
- **效率 = 0.007 / 0.0214 = 32.7%** → 模型只實現了理論上限的 1/3

**結論：不是刪太多 TP，也不是 FP 太少，而是「能觸及的 FP 本來就只有 30.8%」。**

---

## 5. ISM LOH vs LOH.bed 比較

![ISM vs LOH.bed](../../../../../output/synthesis/concluded/loh_verification_analysis/fig3_ism_vs_loh_bed.png)

### 5.1 交叉比較

|  | BED non-LOH | BED LOH | Total |
|--|-------------|---------|-------|
| **ISM non-LOH** | 9,733 | 295 | 10,028 |
| **ISM LOH** | 9,002 | 16,268 | 25,270 |
| Total | 18,735 | 16,563 | 35,298 |

**一致率：73.7%**

### 5.2 不一致的原因

| 類別 | 數量 | 原因 |
|------|------|------|
| **ISM LOH 但 BED non-LOH** (9,002) | 25.5% | ISM 用 HP_Ratio（per-region），LOH.bed 用 phasing block-level 偵測 → ISM 更敏感但也更多雜訊 |
| **BED LOH 但 ISM non-LOH** (295) | 0.8% | 區域內恰好有足夠 minor HP reads → HP_Ratio < 0.9 → ISM 不判為 LOH |

**ISM 的 LOH 定義比 LOH.bed 寬鬆得多**（25,270 vs 16,563 區域）。

### 5.3 AUC 比較

![AUC Comparison](../../../../../output/synthesis/concluded/loh_verification_analysis/fig7_auc_ism_vs_loh_bed.png)

| LOH 定義 | 7-feature LR AUC | 特點 |
|---------|-----------------|------|
| ISM LOH_Strong | **0.912** | 最高，但有 selection bias |
| ISM LOH_Weak | **0.884** | |
| LOH.bed LOH | **0.774** | 覆蓋更廣，AUC 較低 |
| Both LOH (交集) | 0.805 | |
| ISM non-LOH | 0.800 | |
| LOH.bed non-LOH | 0.808 | |

**LOH.bed LOH 的 AUC (0.774) 顯著低於 ISM LOH_Strong (0.912)**，原因：
- LOH.bed 包含大量 ISM 的 LOH_Noise (7,853 區域) 和 LOH_Subclone (1,464)
- 這些子群的 AUC 只有 0.645-0.688
- ISM LOH_Strong 通過 HP_Ratio + VerificationClass 雙重篩選，是高度精選的子群

### 5.4 LOH.bed 的 TP/FP 分布

| 區域 | TP% | FP% | 說明 |
|------|-----|-----|------|
| BED LOH | 83.0% | 17.0% | FP 比例較高 |
| BED non-LOH | 89.3% | 10.7% | FP 比例較低 |

LOH 區域 FP 比例較高（17.0% vs 10.7%），符合預期：LOH 區域的 germline het variants 更容易被誤判為 somatic。

---

## 6. 特徵分布全景

![Feature Distribution](../../../../../output/synthesis/concluded/loh_verification_analysis/fig5_feature_by_loh_subtype.png)

### 6.1 各 LOH_Subtype 的 TP/FP 分離度

| 特徵 | LOH_Strong 分離 | LOH_Weak 分離 | LOH_Noise | None |
|------|----------------|---------------|-----------|------|
| AlleleDelta | TP << FP ✅ | TP < FP ✅ | 微小 | TP < FP |
| HPFineP | TP < FP ✅ | TP < FP | 微小 | 微小 |
| HP_Ratio | TP ≈ FP ❌ | TP ≈ FP | TP ≈ FP | TP ≈ FP |
| CramersV | TP > FP ✅ | TP > FP | TP > FP | TP > FP |
| PairwiseMeanDist | TP ≈ FP | TP ≈ FP | TP ≈ FP | TP ≈ FP |
| LabelHPPermanovaF | 微小 | 微小 | 微小 | 微小 |

**最有區分力的特徵組合：AlleleDelta + CramersV + HPFineP**，且主要在 LOH_Strong/Weak 中有效。

---

## 7. 綜合結論

### 7.1 用戶敘述修正

| 原始敘述 | 修正 |
|---------|------|
| "用 LOH.bed 區域分別" | ISM 自算 HP_Ratio LOH，不讀 LOH.bed |
| "AlleleDelta 是 HP tag 甲基差異" | AlleleDelta 是 ALT-read vs REF-read 甲基化差異 |
| "甲基化新訊號" | 是真實的 ASM 現象，但與 AF confound 高度相關，不超越 AF 本身的資訊量 |
| "可以移除 LOH 下的 germline FP" | 在 LOH_Strong+Weak 子群內可以，但僅覆蓋 30.8% FP |
| "F1 有提升" | 全域 +0.007（微小），但確實是正向的 |

### 7.2 F1 改善小的根因

**主因是覆蓋率 (30.8%)，不是 TP 損失或 FP removal 不足。**

理論上限 ΔF1=+0.021，實際達到 +0.007（效率 33%）。剩餘 69.2% 的 FP 在 LOH_Noise/Subclone/non-LOH 中，這些子群的 AUC 不足以有效過濾。

### 7.3 LOH.bed vs ISM LOH

LOH.bed 更保守（16,563 區域 vs ISM 25,270），但 AUC 更低（0.774 vs 0.805-0.912）。兩者一致率 73.7%。ISM LOH 的優勢來自 HP_Ratio + VerificationClass 的雙重篩選。

---

## 圖片路徑

```
output/synthesis/concluded/loh_verification_analysis/
├── fig1_allele_delta_distribution.png   — AlleleDelta by LOH×TP/FP
├── fig2_hp_ratio_in_loh.png            — HP_Ratio in LOH regions
├── fig3_ism_vs_loh_bed.png             — ISM LOH vs LOH.bed cross-tabulation
├── fig4_f1_breakdown.png               — F1 improvement analysis
├── fig5_feature_by_loh_subtype.png     — Feature distributions by LOH_Subtype
├── fig6_loh_hp_biology.png             — Minor HP reads in LOH (biology)
└── fig7_auc_ism_vs_loh_bed.png         — AUC comparison ISM vs LOH.bed
```
