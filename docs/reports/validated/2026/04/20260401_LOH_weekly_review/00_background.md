<!--
建立時間: 2026-04-01 12:00
更新時間: 2026-04-01 18:00
目標: 提供 LOH 週報審閱的完整研究背景，使讀者從零開始理解 InterSubMod 研究框架、樣本體系、關鍵術語與數據規模
處理範圍:
  - InterSubMod 工具定位與方法概述
  - Region 定義與數據結構
  - LOH-like 操作定義與 LongPhase-TO LOH 的差異
  - 7 個癌症細胞系樣本說明與 cell line 限制
  - Paired / TO 模式定義與差異（含 HP tag 格式詳述）
  - 數據規模與欄位概覽
  - Confound awareness
  - 關鍵術語表
關聯檔案:
  - docs/CURRENT_FOCUS.md
  - docs/experiments/INDEX.md
  - docs/reports/validated/2026/04/20260401_LOH_weekly_review/01_hp_integer_tag_fix.md
  - docs/reports/validated/2026/04/20260401_systematic_observation_O1_O10_summary_01.md
  - docs/reports/validated/2026/03/20260322_TO_FP來源分解摘要_01.md
-->

# 研究背景與數據基礎

**文件定位**：本文件為 LOH 週報審閱系列的第零章，目的是讓任何讀者（包括未接觸過本專案的人）從最基本的概念開始理解整個研究框架。後續章節會假設讀者已閱讀本文件。

---

## 1. InterSubMod (ISM) 是什麼

### 1.1 工具目的

InterSubMod（Inter-Subclonal Methylation Analysis，以下簡稱 ISM）是一個生物資訊分析工具，專門設計用來**從長讀取（long-read）測序數據中偵測腫瘤樣本的亞克隆結構（subclonal structure）**。

要理解這個目的，需要先知道幾個基本背景：

- **腫瘤不是均質的**：一個腫瘤內部往往包含多個不同的細胞群體（稱為「亞克隆」），它們各自帶有不同的突變和表觀遺傳修飾。這種異質性（heterogeneity）是腫瘤演化、抗藥性產生的核心機制。
- **長讀取測序（如 Oxford Nanopore, ONT）**：相比短讀取（Illumina），長讀取可以在**同一條 read 上同時看到**：(a) 體細胞變異（somatic SNV）、(b) 甲基化修飾狀態、(c) 單倍體型（haplotype）歸屬。這三者的共現關係是 ISM 的核心數據基礎。
- **甲基化模式（methylation pattern）**：DNA 上的 CpG 位點可以被甲基化（methylated）或未甲基化（unmethylated）。不同的亞克隆可能有不同的甲基化模式，這種差異稱為表觀遺傳異質性（epigenetic heterogeneity）。

ISM 的做法是：以每一個體細胞變異（somatic SNV）作為錨點（anchor），收集周圍 region 內所有覆蓋該位點的 reads，分析這些 reads 的甲基化模式是否在不同單倍體型之間存在統計顯著的差異。

### 1.2 核心方法

ISM 整合三層資訊來分析每一個 somatic SNV region：

1. **甲基化模式（Methylation Patterns）**：從 BAM 檔案中的 `MM`/`ML` 標籤解碼每條 read 上每個 CpG 位點的甲基化機率。ISM 將甲基化機率 >= 0.5 視為甲基化（methylated），< 0.5 視為未甲基化（unmethylated），並建構一個 read-by-CpG 的二元矩陣。
2. **體細胞變異（Somatic SNVs）**：從 VCF 檔案讀取 somatic caller（如 ClairS 或 ClairS-TO）偵測到的變異位點，每個位點作為一個分析 region 的中心。ISM 判斷每條 read 在該位點支持 REF allele 還是 ALT allele。
3. **單倍體型（Haplotypes）**：從 BAM 檔案的 `HP` 標籤判斷每條 read 屬於哪個 haplotype（HP:1 或 HP:2）。這個標籤由上游的 phasing 工具（LongPhase）寫入。

ISM 使用多種距離度量（L1 / L2 / NHD / Bernoulli / Jaccard）計算 reads 之間的甲基化距離，並透過統計檢驗（Fisher's Exact Test / Chi-square / PERMANOVA）評估分群結果的顯著性。最後，ISM 對每個 region 產出一組特徵值（如 `VerificationClass`、`Quality_Score`、`PairwiseMedianDist`、`CramersV`、`AlleleDelta` 等），作為該 somatic SNV 的 read-level epigenetic evidence。

### 1.3 在 variant calling pipeline 中的定位

**ISM 不是一個 variant caller**。它不負責「呼叫」（call）哪些位置有 somatic mutation。那是上游工具（如 ClairS、ClairS-TO、DeepSomatic）的工作。

ISM 的角色是**提供額外的 read-level epigenetic context**，幫助下游決策判斷一個已經被 caller 呼叫出來的 variant 是否可信（true positive, TP）還是假陽性（false positive, FP）。換句話說：

- Caller 告訴你「這個位置可能有突變」
- ISM 告訴你「覆蓋這個位置的 reads，它們的甲基化模式和單倍體型分佈，是否支持這是一個真實的腫瘤突變」

這個定位非常重要，因為它決定了 ISM 的研究方向：**我們不是在發明新的 calling 演算法，而是在挖掘 read-level 的表觀遺傳學證據，看它能否提升現有 caller 的精確度。**

---

## 2. Region 定義

### 2.1 什麼是一個 Region

ISM 中的**每一個 Region 對應一個 somatic variant**（1:1 關係）。具體而言：

- 以 somatic caller 報告的每一個 SNV 的基因組座標為中心
- 向兩側各延伸 `window_size` bp（預設 ±500 bp），形成一個分析窗口
- ISM 收集該窗口內所有 reads 的甲基化、allele support、haplotype 資訊
- 對這些 reads 進行聚合統計，產出**一行**（row）完整的特徵向量

因此，master dataset 中的**每一行 = 一個 Region = 一個 somatic variant = 一個 TP 或 FP**（依據 truth set 比對判定）。

### 2.2 Region 內的 read 聚合

Region 層級的統計是對該窗口內所有 reads 的聚合結果。例如：

- `NumReads`：窗口內的總 read 數
- `HP1FamilyN`：屬於 HP1 家族（HP1 + HP1-1）的 read 數
- `HP_Ratio`：HP1 家族 reads 佔有效 HP reads 的比例
- `region_methyl_mean`：窗口內所有 CpG 的平均甲基化水平

這意味著 Region 層級的特徵是**近似值**——它反映的是窗口範圍內 reads 的統計摘要，而非單一位點的精確狀態。

---

## 3. LOH-like 定義

### 3.1 ISM 的 LOH-like 操作定義

**LOH（Loss of Heterozygosity，雜合性喪失）** 是指原本在兩條染色體上各有不同 allele 的位置，因為其中一條染色體的片段丟失或被另一條複製取代，導致該區域只剩下一個 allele 的狀態。

ISM 使用 `HP_Ratio` 作為 region 層級的 LOH 近似指標（proxy）：

- **LOH-like 判定條件**：`HP_Ratio < 0.1` 或 `HP_Ratio > 0.9`

**HP_Ratio 計算公式**（來源：`src/core/RegionProcessor.cpp`）：

```
HP_Ratio = (HP1FamilyN + ε) / (HP1FamilyN + HP2FamilyN + 2ε)
```

其中 ε = 0.001（避免除以零）。

- `HP1FamilyN` = HP1 家族的 read 數（包含 HP1 + HP1-1）
- `HP2FamilyN` = HP2 家族的 read 數（包含 HP2 + HP2-1）
- 不包含 HP0（未分配）和 HP3（ambiguous somatic）

當 `HP_Ratio` 接近 0.5 時，兩個 haplotype 大致平衡（non-LOH）；當極端偏向一側時，標記為 LOH-like。

### 3.2 與 LongPhase-TO LOH 的差異

ISM 的 LOH-like 與上游 LongPhase-TO 報告的 LOH 是**不同層級的定義**，兩者不完全一致：

| 面向 | ISM LOH-like | LongPhase-TO LOH |
|:---|:---|:---|
| 分析層級 | Region-level（單一 variant 窗口內 reads 聚合） | Phase block-level（基因型推論） |
| 資訊來源 | HP tag read 計數比例 | LOH.bed + phased VCF 基因型 |
| 定義方式 | HP_Ratio < 0.1 or > 0.9（經驗閾值） | 基因型推論顯示單倍體型丟失 |
| 精確度 | 近似值（受窗口大小、read sampling 影響） | 較精確（基於整個 phase block） |

### 3.3 原子式分析的思考

目前的分析使用 **region-level LOH-like** 作為整體分析單元。一個更精確的替代方案是**逐位點判斷**：對每個 variant 位點，直接查詢 LongPhase-TO 的 LOH.bed 和 phased VCF，判斷該位點是否落在 LOH block 內。這種原子式分析更接近 LongPhase-TO 的原始定義，但目前尚未實作。

Region-level LOH-like 的優勢在於不依賴外部 LOH annotation，可以純粹從 ISM 自身的 read 資訊推導；缺點是作為近似值，會產生假陽性（非 LOH 但 HP imbalance）和假陰性（LOH 但窗口內仍有足夠雙側 reads）。

---

## 4. 本週研究核心問題

本週（2026-03-25 至 2026-03-31）的研究聚焦於兩個核心問題：

### 4.1 LOH 相關特徵能否用於鑑別 FP？

核心問題是：LOH-like 的 region 中，FP 的比例是否特別高？如果是，那麼「LOH-like」就可以作為過濾 FP 的特徵。

### 4.2 Paired vs TO 的本質差異

在深入 LOH 分析之前，必須先理解 Paired 和 TO 是**本質上不同的問題空間**：

1. **候選集（universe）差異**：Paired 模式有 normal 對照，FP rate 僅 ~1.04%；TO 沒有 normal，FP rate 高達 ~30.6%。兩者的候選集組成從根本上不同。
2. **FP 組成差異**：Paired FP 是 caller 技術性誤報；TO FP 的 ~99% 是 germline variant 被錯認為 somatic（詳見第 6.2 節）。
3. **Phasing 品質差異**：Paired 使用 LongPhase-S（高品質 germline SNP phasing）；TO 使用 LongPhase-TO（依賴腫瘤自身 variant，品質較差）。
4. **特徵方向差異**：LOH enrichment 方向完全相反——Paired 中 LOH-like 傾向 FP 富集（enrichment > 1）；TO 中 LOH-like 反而傾向 TP 富集（enrichment < 1）。

這些差異意味著 **Paired 和 TO 必須分離處理**——不能混合訓練模型或使用相同的 feature set，否則 TO 的高 FP rate 和反向特徵方向會汙染 Paired 的結果。這是本週多份報告的核心前提。

---

## 5. 七個樣本說明

### 5.1 Cell Line 與臨床樣本的差異

本研究使用的 **7 個樣本全部是癌症細胞系（cell line）**，不是臨床患者組織樣本。每個樣本都同時以 Paired 和 TO 兩種模式運行，產生 14 個 dataset。

**為什麼使用 cell line：**

- Cell line 是方法開發和基準測試的**標準做法**。SEQC2 consortium、DeepSomatic、Lancet2 等主流工具和計畫都使用 cell line 建立 truth set 和評估效能。
- Cell line 具有穩定可重複的特性，能提供 SEQC2 金標準 truth set，使 TP/FP 判定有明確依據。

**Cell line 的已知限制：**

- **缺乏腫瘤微環境**：cell line 經長期體外培養，缺少臨床腫瘤中的免疫細胞浸潤、基質細胞、血管結構等微環境成分。
- **Purity 接近 100%**：cell line 幾乎是純腫瘤細胞，而臨床樣本通常混合正常細胞（purity 20%-80%），低 purity 會影響 allele frequency 和 phasing 品質。
- **長期培養偏差**：cell line 可能因長期培養累積額外突變或甲基化漂移，偏離原始腫瘤特徵。
- **甲基化特異性**：cell line 的全基因組甲基化水平通常低於原始組織，可能影響 ISM 甲基化分析的外部效度。

**推廣到臨床時的需求**：tissue-origin correction、purity-aware analysis、matched normal availability 評估。這些校正在本研究範圍之外，但對未來臨床應用是必要的。

### 5.2 樣本總覽

| 樣本名稱 | 癌症類型 | 基因組特徵 | 用途角色 | 說明 |
|:---|:---|:---|:---|:---|
| HCC1395 5kHz | Triple-negative breast cancer (TNBC) | BRCA1/TP53 突變，高度基因組不穩定 | Discovery（主要開發樣本） | 使用 5kHz sampling rate 的 ONT 定序；是整個專案最早、最深入分析的樣本，大部分方法論都在此樣本上開發和驗證 |
| HCC1395_DORADO | Triple-negative breast cancer (TNBC) | 同 HCC1395，不同 basecalling 平台 | Platform validation | 使用 DORADO basecaller 重新處理的同一樣本；用於驗證 ISM 結果是否因 basecalling 平台不同而改變 |
| COLO829 | Melanoma（黑色素瘤） | BRAF V600E 突變，UV 相關突變譜 | Cross-cancer validation | SEQC2 consortium 的標準參考樣本，有公開的 truth set；突變數量中等 |
| H1437 | Lung adenocarcinoma（肺腺癌） | KRAS 突變 | Cross-cancer validation | 肺癌代表樣本之一；TO 模式下 FP rate 較低 |
| H2009 | Lung adenocarcinoma（肺腺癌） | KRAS/TP53 突變 | Heterogeneity study | TO 模式下 FP rate 極高（74.6%），是其他樣本的 2-5 倍。佔整體 TO FP 的 36.2%，是異質性研究的重要對象 |
| HCC1937 | BRCA1-mutant breast cancer | BRCA1 germline 突變，HRD 表型 | Cross-cancer validation | DNA repair deficiency 代表樣本 |
| HCC1954 | HER2+ breast cancer | HER2 擴增，TP53 突變 | Cross-cancer validation | 高度基因組重排，LOH 區域特別多 |

### 5.3 圖示：各樣本/模式的 Region 數量

![各樣本/模式的 region 數量](../figures/20260401_observation_figures/O01/fig03_sample_mode_region_counts.png)

### 5.4 圖示：TP/FP 比例

![TP/FP 比例](../figures/20260401_observation_figures/O01/fig04_truth_label_composition.png)

---

## 6. Paired vs TO 模式定義

### 6.1 Paired 模式

**定義**：同時擁有腫瘤（tumor）和配對正常樣本（matched normal）的分析模式。

- **Somatic calling**：使用 ClairS（paired somatic caller），透過腫瘤與正常樣本的比較來識別體細胞變異。因為有正常樣本作為對照，能有效過濾掉 germline variant，false positive rate 很低。
- **Phasing**：使用 LongPhase-S（standard phasing），利用 germline heterozygous SNP 進行 read phasing。Phasing 品質高，HP tag 分配準確。
- **HP tag 格式**：字串格式 `HP:Z:`，包含以下值：
  - `HP:Z:1` — germline haplotype 1
  - `HP:Z:2` — germline haplotype 2
  - `HP:Z:1-1` — somatic-supported haplotype 1（read 同時攜帶 somatic ALT allele 且被追溯歸屬到 germline HP1）
  - `HP:Z:2-1` — somatic-supported haplotype 2（同上，歸屬到 HP2）
  - `HP:Z:3` — ambiguous somatic（read 攜帶 somatic ALT 但無法唯一歸屬到 HP1 或 HP2）
  - 無 HP tag — unphased reads（HP0），phasing 工具無法將該 read 歸入任何 haplotype
- **說明**：「1-1」中的第二個數字代表 somatic support 標記，**不是** phase block ID。`HP:Z:1-1` 意為「屬於 HP1 家族，且有 somatic allele 支持」。
- **FP rate**：平均約 **1.04%**（每 100 個被呼叫的 variant，約 1 個是假陽性）。

### 6.2 Tumor-Only (TO) 模式

**定義**：僅有腫瘤樣本，沒有配對正常樣本的分析模式。

- **Somatic calling**：使用 ClairS-TO（tumor-only somatic caller）。因為沒有正常樣本對照，caller 必須依靠 population database（如 gnomAD、PON）來過濾 germline variant，但仍會有大量 germline variant 被錯誤地呼叫為 somatic（造成高 FP rate）。
- **Phasing**：使用 LongPhase-TO（tumor-only phasing），利用腫瘤樣本自身的 heterozygous variant 進行 phasing。品質低於 Paired 模式。
- **HP tag 格式**：整數格式 `HP:i:`，包含以下值：
  - `HP:i:1` — germline HP1（等同 Paired 的 `HP:Z:1`）
  - `HP:i:2` — germline HP2（等同 `HP:Z:2`）
  - `HP:i:11` — somatic-supported HP1（等同 `HP:Z:1-1`）
  - `HP:i:21` — somatic-supported HP2（等同 `HP:Z:2-1`）
  - `HP:i:33` — ambiguous somatic（等同 `HP:Z:3`，**不是第三個 haplotype**）
  - 無 HP tag — unphased reads（HP0）
- **ISM 內部統一處理**：ISM 的 ReadParser（`src/core/ReadParser.cpp`）在讀取 BAM 時，將 TO 整數格式自動轉換為與 Paired 相同的字串格式（`11 -> "1-1"`, `21 -> "2-1"`, `33 -> "3"`），下游分析使用統一介面。這也是 HP Integer Tag Bug 的根源——修正前 ISM 未正確辨識整數格式的 HP tag（詳見 `01_hp_integer_tag_fix.md`）。
- **FP rate**：平均約 **30.6%**（每 100 個被呼叫的 variant，約 31 個是假陽性），是 Paired 模式的 **~30 倍**。

### 6.3 TO FP 的本質

TO 的 30.6% FP 中，根據 provenance 分析：

- **~98.55%-98.74% 是 `paired_raw_absent`**：這些 variant 在 paired 模式下根本不會被呼叫為 somatic——它們絕大多數是 germline variant，因為缺少 normal 樣本對照而無法被過濾。
- **<1% 是 caller 技術性誤報**（如 `paired_persistent_final_fp`，HCC1395=87 個，HCC1395_DORADO=77 個）。

換言之，TO 的高 FP rate 本質上是**「缺失 normal 導致 caller candidate universe 擴張」**的問題——caller 的搜尋空間變大了，而非 caller 本身品質下降。這使得 ISM 甲基化特徵對 TO FP 過濾的效用極為有限，因為 germline variant 本身就會展現正常的 allele-specific methylation（ASM）模式，與 somatic TP 的甲基化特徵難以區分。

### 6.4 為什麼分離處理是必要的

綜合以上，Paired 和 TO 的區分貫穿本研究的每一個環節：

1. **FP 組成完全不同**：Paired 的 FP 主要是 caller 的技術性誤報；TO 的 FP 大多是 germline variant。這兩類 FP 在甲基化特徵上表現不同。
2. **LOH 方向相反**：Paired LOH 傾向 FP 富集，TO LOH 傾向 TP 富集（第 4.2 節）。
3. **Phasing 品質不同**：TO phasing 的 HP assignment 受限於 somatic allele 造成的 HP imbalance，導致 LOH-like 的成因與 Paired 不同。
4. **特徵有效性不同**：系統性觀察（O1-O10）反覆確認，在 TO 模式下幾乎沒有單一 ISM 特徵能有效區分 TP/FP（所有 AUC < 0.58），但在 Paired 模式下 `GQ`（caller quality）AUC 可達 0.811。
5. **研究策略必須分離**：Paired 和 TO 不能混合訓練模型或使用相同的 feature set，否則 TO 的高 FP rate 和反向特徵方向會汙染 Paired 的結果。

---

## 7. 數據規模

本研究使用的主資料集（master dataset）基本規格如下：

| 項目 | 數值 |
|:---|:---|
| 總行數（rows） | 748,391 |
| 總欄位數（columns） | 116 |
| 資料來源 | 7 samples x 2 modes = 14 datasets 的 `significance_summary.csv` 合併 |
| 資料狀態 | Post HP-fix（2026-03-30 之後，已修正 TO HP integer tag bug） |
| 儲存格式 | `all_region_rows.tsv.gz`（壓縮 TSV） |
| 每行代表 | 一個 somatic SNV region 的完整特徵向量（一行 = 一個 TP 或 FP） |
| 主要欄位類別 | 基本資訊（chr, pos, ref, alt）、caller 特徵（AF, GQ, DP）、ISM 特徵（VerificationClass, QS, PairwiseMedianDist, CramersV, AlleleDelta 等）、HP/LOH 特徵（HP_Ratio, effective_hp_reads, Potential_LOH 等）、truth label（TP/FP） |

---

## 8. Confound Awareness

### 8.1 已確認的 confound 問題

系統性觀察（O11, O12, O13）反覆確認，許多看似有效的 ISM 特徵實際上由 confound 驅動：

- **n_reads（NumReads）是最大的 confound**：TP 平均有更多 reads（paired 中差距尤其明顯），導致依賴 read count 的衍生特徵（如 epipolymorphism、effective_hp_reads）呈現虛假的區分力。
  - 典型案例：epipolymorphism AUC = 0.845 → 殘差化 n_reads 後 AUC = 0.530（近隨機）
- **AF（allele frequency）是 LOH 分析的 confound**：AF 與 HP_Ratio 高度相關，導致 LOH 特徵的表現部分來自 AF 的 confound 而非 LOH 本身（O12 已確認 AlleleDelta = AF confound）。
- **L2 collider bias**：將近常數特徵 residualize on AF 時，可能產生虛假信號。必須以 L3（AF-bin 交叉驗證）確認（O12 發現）。

### 8.2 閱讀 AUC 數字的注意事項

閱讀本系列報告中任何 AUC 數字時，應確認：

1. 該 AUC 是否已控制主要 confound（n_reads, AF）？
2. 是否區分 Paired / TO？（混合模式的 AUC 無意義）
3. 不同樣本的 confound 程度不同——全樣本統計可能被某些特殊樣本（如 H2009）主導。

### 8.3 目前缺乏的示意圖

以下示意圖有助理解但目前尚未製作：

- LOH 位點的 reads HP 分佈示意圖（正常 vs LOH-like 對比）
- Phasing 原理示意圖（LongPhase-S vs LongPhase-TO 的 HP assignment 流程）

---

## 9. 關鍵術語表

以下術語在整個 LOH 週報系列中頻繁出現，請務必理解其定義：

### 9.1 HP tag (Haplotype tag)

BAM 檔案中記錄每條 read 屬於哪個 haplotype 的標籤。由上游的 phasing 工具（LongPhase）寫入。

**Paired 模式（LongPhase-S）**——字串格式 `HP:Z:`：

| Tag 值 | 語義 | ISM 家族歸屬 |
|:---|:---|:---|
| `HP:Z:1` | Germline haplotype 1 | HP1-family |
| `HP:Z:2` | Germline haplotype 2 | HP2-family |
| `HP:Z:1-1` | Somatic-supported HP1（read 攜帶 somatic ALT 且追溯到 germline HP1） | HP1-family |
| `HP:Z:2-1` | Somatic-supported HP2 | HP2-family |
| `HP:Z:3` | Ambiguous somatic（無法唯一歸屬到 HP1 或 HP2） | 排除（NHP3 計數） |
| 無 HP tag | Unphased（phasing 工具無法歸類） | 排除（NHP0 計數） |

**TO 模式（LongPhase-TO）**——整數格式 `HP:i:`：

| Tag 值 | 等同 Paired | 語義 |
|:---|:---|:---|
| `HP:i:1` | `HP:Z:1` | Germline HP1 |
| `HP:i:2` | `HP:Z:2` | Germline HP2 |
| `HP:i:11` | `HP:Z:1-1` | Somatic-supported HP1 |
| `HP:i:21` | `HP:Z:2-1` | Somatic-supported HP2 |
| `HP:i:33` | `HP:Z:3` | Ambiguous somatic（**不是**第三個 haplotype） |
| 無 HP tag | 無 HP tag | Unphased (HP0) |

**重要澄清**：
- `HP:Z:1-1` 中的 `-1` 代表 **somatic support 標記**，不是 phase block ID。
- `HP:i:33` / `HP:Z:3` 不代表第三個 haplotype，而是「read 攜帶 somatic allele 但無法確定歸屬哪個 haplotype」的 ambiguous 狀態。

### 9.2 HP 家族（HP Family）

ISM 將 reads 依 HP tag 分為兩個家族進行分析：

- **HP1-family**：包含 `HP:Z:1`（germline HP1）+ `HP:Z:1-1`（somatic-supported HP1）
- **HP2-family**：包含 `HP:Z:2`（germline HP2）+ `HP:Z:2-1`（somatic-supported HP2）
- **排除**：`HP:Z:3`（ambiguous somatic）、`HP:Z:0` / 無 tag（unphased）

HP1-1 和 HP2-1 的存在使 ISM 能區分「純 germline phasing 支持的 reads」和「同時有 somatic allele 證據的 reads」，但在 merged 分析（Stage 1）中，同家族的 reads 被合併計算。

### 9.3 HP0 (Unphased Reads)

沒有 HP tag 的 reads，表示 phasing 工具無法將該 read 歸入任何 haplotype。在 ISM 的 HP 分析中被排除，不計入 effective_hp_reads。HP0 比例高通常代表該 region 的 phasing 品質差。

### 9.4 HP3 (Ambiguous Somatic Reads)

HP tag 為 `HP:Z:3` 或 `HP:i:33` 的 reads。這些 reads 攜帶 somatic allele support，但 phasing 工具無法唯一確定它們屬於 HP1 還是 HP2。

**重要**：HP3 **不是第三個 haplotype**，而是「ambiguous somatic」的標記。在 ISM 的 HP 分析中被排除，不計入 effective_hp_reads。HP3 也不是 LOH 的標記。

### 9.5 effective_hp_reads (eff_hp)

一個 region 中屬於 HP1-family 和 HP2-family 的 reads 總數：

```
effective_hp_reads = HP1FamilyN + HP2FamilyN
```

- **包含**：HP1 + HP1-1 + HP2 + HP2-1（兩個 haplotype 家族的所有 reads）
- **排除**：HP0（unphased，無 HP tag）和 HP3（ambiguous somatic，無法確定 haplotype）

語義：衡量該 region 有多少 reads 可以被用於 haplotype-level 分析。

- `eff_hp = 0` → 該 region 完全沒有 phasing 資訊，無法判斷 LOH 狀態
- `eff_hp` 越高 → HP-based LOH 信號越可靠

**重要分佈差異**：
- Paired：TP median eff_hp = 68, FP median = 32（差距 2.1 倍，AUC = 0.727）
- TO：TP median = 72, FP median = 63（差距僅 1.1 倍，AUC = 0.564）

### 9.6 LOH-like

一個 region 中 reads 的 haplotype 分佈極度不平衡的狀態。操作型定義：

- `HP_Ratio < 0.1` 或 `HP_Ratio > 0.9`

LOH-like 不一定是真正的 LOH（可能是 phasing 造成的 HP imbalance），但它標記了在 haplotype 層面上資訊量不足的 region。

### 9.7 HP_Ratio 與 LOH-like 的關係

`HP_Ratio` 是連續值（0 到 1），`LOH-like` 是根據 HP_Ratio 的二元判定。HP_Ratio 反映的是 region 層級的 HP 平衡度，而非基因型層級的真實 LOH 狀態（參見第 3.2 節）。

### 9.8 region_methyl_mean

Region 內所有 CpG 位點的平均甲基化水平（0 到 1）。反映該窗口的整體甲基化狀態——高值（接近 1）表示高度甲基化區域，低值（接近 0）表示低甲基化區域。

### 9.9 Tier 分層

根據 `effective_hp_reads` 將 regions 分為不同信賴層級：

| Tier | effective_hp_reads 範圍 | 語義 |
|:---|:---|:---|
| C0 | = 0 | 完全無 HP 資訊，LOH 判斷不可能 |
| C | 1 - 9 | HP 資訊不足，LOH 判斷極不穩定 |
| B | 10 - 29 | HP 資訊中等，LOH 判斷有參考價值但需謹慎 |
| A | 30 - 49 | HP 資訊充足，LOH 判斷可信 |
| A+ | >= 50 | HP 資訊非常充足，LOH 判斷高度可信 |

**關於 Tier 邊界**：

- 這些邊界是**經驗性的**（非理論推導），基於觀察 enrichment 在不同 eff_hp 水平下的穩定性選定。
- 實際分佈數據支持這些邊界具有區分力：
  - Paired 中，**TP 的 67.62% 落在 A+ tier**，而 FP 只有 25.66%——差距顯著。
  - TO 中，TP 和 FP 在各 tier 的分佈差距較小，反映 TO 模式下 eff_hp 的區分力下降。
- Tier 的主要用途：區分「有足夠 phasing support 的可信 LOH 信號」vs「read 數不足導致的 noise」。
- 未來可根據 read 深度分佈更精細調整，例如引入 sample-specific 或 coverage-normalized 的 tier 邊界。

HP-fix 之前，TO 模式有大量 regions 落在 C0/C tier（因為 HP tag 未被正確解析）；修正後，88% 的 TO regions 進入 Tier A/A+。

### 9.10 Enrichment

Enrichment 是衡量某個特徵（如 LOH-like）在 FP 與 TP 之間偏好程度的指標：

```
enrichment = FP_LOH% / TP_LOH%
```

- **enrichment > 1**：LOH-like 中 FP 的比例高於 TP → FP 富集 → LOH 可作為 FP 風險標記
- **enrichment < 1**：LOH-like 中 TP 的比例高於 FP → TP 富集 → LOH 不能用來過濾 FP
- **enrichment = 1**：LOH-like 對 FP/TP 沒有區分力

### 9.11 VerificationClass

ISM 根據統計檢驗結果和 haplotype/allele 分佈模式，對每個 region 給出的分類標籤：

| 類別 | 語義 |
|:---|:---|
| Strong | 強證據支持該 variant 是真實的 subclonal methylation difference |
| Weak | 弱證據，有甲基化差異但統計顯著性不足 |
| Subclone | 偵測到 subclonal 結構的證據 |
| Noise | 沒有有意義的甲基化差異，可能是噪聲 |

### 9.12 QS (Quality Score)

ISM 對每個 region 計算的綜合品質評分，結合了多個特徵（統計顯著性、效應量、HP 一致性等）。QS 的設計目的是提供一個單一數值來判斷 variant 的可信度。

**重要警示**：系統性觀察已確認，TO 模式下 QS 完全失效（AUC = 0.497，等同隨機），主要原因是 LOH penalty 在 TO 中方向反轉（LOH 在 TO 是 TP enrichment，但 QS 把它當作負面指標）。因此，TO 的 QS 需要根本性重設計（已在 `07_qs_mode_aware_change.md` 中實施第一步修正）。

### 9.13 AUC (Area Under ROC Curve)

衡量一個特徵區分 TP 和 FP 能力的標準指標：

- **AUC = 0.5**：等同隨機猜測，該特徵完全沒有區分力
- **AUC = 0.6-0.7**：弱區分力
- **AUC = 0.7-0.8**：中等區分力
- **AUC > 0.8**：強區分力
- **AUC = 1.0**：完美區分

在本研究中，我們關注的是各個 ISM 特徵對 TP/FP 的 AUC。一個反覆確認的結論是：TO 模式下所有 ISM 特徵的 AUC < 0.58（接近隨機），而 Paired 模式下某些 caller 特徵（如 GQ）可達 0.811。**閱讀任何 AUC 數字時，請參照第 8 節的 confound awareness 注意事項。**

---

## 待驗證問題（已驗證 / 已更新）

### 已解決

1. **~~LOH-like 的成因分離~~** ✅ C13 確認：TO LOH-like 主要由 phasing HP imbalance 造成（非真實基因組 LOH）。39,724 個 Paired→TO 翻轉位點中 71.6% 的 TO min(HP1,HP2)=0。TO 的「LOH-like」本質上偵測 somatic clonality。

4. **~~Tier enrichment 的穩定性~~** ✅ C12 確認：6/7 樣本 A+ 實際優於 A（Simpson's Paradox）。Per-sample 驗證已完成，非少數樣本驅動。

### 尚未解決

2. **跨樣本穩定性（H2009 主導效應）**：H2009 TO FP rate 74.6%。需 leave-one-out 驗證。定位 P2。
3. **FN 的 LOH 分佈**：需 FN ISM 數據。對應 P2「O9 FN 觀察」。
5. **Cell line 到臨床的可遷移性**：purity ~100% vs 20-80%。定位 P2。
6. **Confound 控制後的殘餘信號**：需 DAG 分析。定位 P2。

---

## 認知門檻補充建議

以下是審閱後續章節時可能遇到的認知門檻，建議讀者預先理解：

1. **Phasing 的局限性**：Phasing 不是完美的。即使 HP tag 被正確寫入，phasing 本身的錯誤（如 switch error）也會導致 reads 被錯誤分配到另一個 haplotype。這種錯誤在 TO 模式中更常見。
2. **Enrichment 不等於因果**：enrichment > 1 表示 FP 在 LOH-like 中比例偏高，但不代表 LOH 導致 FP。可能是第三個變數（如 low coverage、特定基因組區域）同時造成 LOH-like 和 FP 偏高。
3. **Region-level proxy 的限制**：本報告中的 LOH-like 是 region 層級的 HP_Ratio 近似值，不等於真正的基因型 LOH。解讀時需注意這一層級差異（參見第 3 節）。
4. **Multi-testing 問題**：748,391 個 regions 的獨立檢驗會產生大量 false discovery。目前尚未對所有統計檢驗進行系統性的多重比較校正（FDR）。
5. **Confound 無處不在**：系統性觀察（O11, O12）已確認，許多看似有效的特徵（如 epipolymorphism AUC=0.845）實際上是 confound（如 n_reads）造成的虛假信號。殘差化（residualization）或分層（stratification）後 AUC 大幅下降至接近隨機。閱讀任何 AUC 數字時，都應確認是否已控制 confound（參見第 8 節）。
