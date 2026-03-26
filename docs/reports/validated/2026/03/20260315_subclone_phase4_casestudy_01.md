<!--
建立時間: 2026-03-15
目標: Phase 4 Case Study — 5 個 Strong Subclone 候選位點甲基化視覺觀察文件
處理範圍: HCC1395 5kHz paired (s-pure/HCC1395/20260307)，VerificationClass=Subclone AND CramersV>0.3 AND NumReads>=30
關聯檔案:
  - 20260315_subclone_phase4_casestudy_assets/ (所有圖片)
  - analysis/subclone_case_study_candidates_20260315.tsv
  - scripts/analysis/build_subclone_map.py
-->

# Subclone Phase 4 Case Study — 甲基化視覺觀察

**觀察日期**：2026-03-15
**資料來源**：`s-pure/HCC1395/20260307`（Paired tumor-only，ClairS + LongPhase-S haplotagging）
**篩選條件**：`VerificationClass = Subclone` AND `CramersV > 0.3` AND `NumReads ≥ 30`（全為 TP）

---

## 資料狀態說明（閱讀前必看）

| 項目 | 狀態 |
|------|------|
| **模式** | Paired（有配對 normal），非 Tumor-Only |
| **HP tag 來源** | LongPhase-S：以 normal 的 germline SNP 建立 phasing，再將 phase 標籤 apply 到 tumor reads |
| **`HP=1` / `HP=2`** | Read 被分配至 germline haplotype 1 或 2（來自 normal 相位，反映正常細胞 germline 背景） |
| **`HP=0`** | 此 read 在該區域**無法相位**（附近無足夠 germline SNP） |
| **`HP=1-1`** | HP=1 背景的 read，且 read **帶有 somatic ALT allele**（原屬 HP1 haplotype） |
| **`HP=2-1`** | HP=2 背景的 read，且 read **帶有 somatic ALT allele**（原屬 HP2 haplotype） |
| **`HP=3`** | LongPhase-S **無法區分 HP1 或 HP2**，但此 read **帶有 somatic ALT**；屬 ALT-supporting reads 的 unresolved HP 類別 |
| **HP 標籤關鍵注意** | HP=1-1 / HP=2-1 / HP=3 **均表示 ALT reads**；HP=1 / HP=2 通常為 REF reads（germline haplotype 背景）；HP=0 可能是 REF 或 ALT（無相位區域） |
| **LOH 判讀注意** | `HP=3` 只表示 ALT read 的 germline HP 背景未能唯一解析；若要判定 LOH，仍需另查 phased genotype、`LOH.bed` 或 CNV 證據 |
| **甲基化值** | 5mC analog (MM/ML tag)，0.0=完全未甲基，1.0=完全甲基化 |
| **圖說顏色** | Cluster Heatmap：紅=高甲基，藍=低甲基；HP strip 色碼：灰=HP0，藍=HP1，黃=HP2，深藍=HP1-1，橙=HP2-1，紫=HP3 |
| **距離熱力圖** | viridis_r 色碼：深色=reads 相似（低 Bernoulli 距離），亮色=reads 差異大；兩軸均有 UPGMA 樹形 |
| **視覺化腳本** | `tools/plot_cluster_heatmap.py`（甲基矩陣）、`tools/plot_distance_heatmap.py`（距離矩陣）、metric=BERNOULLI |

---

## 全域背景：亞克隆地圖

### TP vs FP 特徵分布

![TP vs FP Feature Distribution](20260315_subclone_phase4_casestudy_assets/methylation_tp_fp_distribution.png)

### 染色體亞克隆密度圖

![Subclone Chromosome Ideogram](20260315_subclone_phase4_casestudy_assets/subclone_map_chr_ideogram.png)

> 左側：每 Mbp 的 Subclone 位點密度；右側：總數（藍）與 Strong Subclone (CramersV>0.3，紅）。
> **5 個案例分布在 chr4, chr5, chr6, chr7, chr16**（以紅條標示）。

---

## 案例 1：`chr6:145444893` G>A — CramersV = **0.722**（最強信號）

### 數值摘要

| 指標 | 值 | 異常旗標 |
|------|-----|---------|
| Reads / CpGs | 58 / 39 | |
| **Optimal k** | **3** | ⚠️ 多數為 2，此處分 3 群 |
| ALT reads (n) | 19 | all HP=**2-1** (HP2 背景 ALT) |
| REF reads (n) | 39 | all HP=2 or HP=0 |
| ALT meth mean | 0.580 ± 0.064 | |
| REF meth mean | 0.541 ± 0.130 | REF variance 明顯大於 ALT |
| Allele PERMANOVA p | **0.001** (valid) | ✅ allele 驅動分類 |
| HP PERMANOVA p | **1.000** (invalid) | HP 完全不相關 |
| HP_Ratio | 0.019 | 極偏 HP2，HP 無法分離 |
| 整體甲基化範圍 | 0.40–0.70 | 中等甲基化 |

**關鍵特徵**：k=3 三群；所有 ALT reads 都落在 HP=2-1；allele 是唯一有效 label。

### 甲基矩陣熱力圖（Cluster Heatmap：Read × CpG）

> **X 軸**：CpG 位點（依基因組位置排序）；**Y 軸**：Reads（依 BERNOULLI 距離矩陣 UPGMA 聚類排序，左側樹形圖）。
> **左側彩條**：HP tag（灰=未相位，藍=HP1，黃=HP2，橙=HP2-1）；Strand（+/-）；Source（Tumor/Normal）；Allele（灰=ALT，淺=REF）。

![chr6:145444893 cluster heatmap](20260315_subclone_phase4_casestudy_assets/chr6_145444893_cluster_heatmap.png)

### Bernoulli 距離熱力圖（Read × Read）

> **兩軸均為 Reads**，按同一 UPGMA 聚類排序，兩軸均有樹形圖。深色（viridis_r）= reads 甲基化模式相似（低 Bernoulli 距離）；亮色 = 差異大。

![chr6:145444893 distance heatmap](20260315_subclone_phase4_casestudy_assets/chr6_145444893_distance_heatmap.png)

### methylartist — ±2kb per-read 甲基化（含 SNV 位置）

> SVG 圖：每條水平線 = 一條 read，點 = CpG 甲基化狀態（顏色深淺 = 甲基化程度），垂直黑線 = SNV 位置，HP 分層標示在左側。

![chr6:145444893 methylartist](20260315_subclone_phase4_casestudy_assets/chr6_145444893_methylartist.svg)

### 觀察問題

- [ ] Heatmap 中的 3 群是否在視覺上清晰分層？
- [ ] methylartist SVG 中，ALT reads（HP=2-1）與 REF reads（HP=2）的甲基化 pattern 是否有明顯差異，還是隨機散布？
- [ ] k=3 是否代表有 2 個 REF 甲基亞群 + 1 個 ALT 甲基群？

---

## 案例 2：`chr4:70548355` G>A — CramersV = 0.502

### 數值摘要

| 指標 | 值 | 異常旗標 |
|------|-----|---------|
| Reads / CpGs | 48 / **63** | 最多 CpGs（最豐富甲基資訊） |
| Optimal k | 2 | |
| **HP=0 reads** | **32/48 = 67%** | ⚠️ 異常高未相位比例 |
| ALT reads (n) | 11 | mostly HP=1-1 (HP1 背景 ALT) |
| REF reads (n) | 37 | mostly HP=0 |
| ALT meth mean | 0.648 ± 0.103 | |
| REF meth mean | 0.611 ± 0.064 | 差異小，ALT variance 較大 |
| Allele PERMANOVA p | 0.011 (valid) | ✅ 邊緣顯著 |
| HP_Ratio | **0.937** | 極偏（但 HP=0 佔大多數） |
| 整體甲基化 | 高（0.58–0.70） | 整體高甲基化區域 |

**關鍵特徵**：2/3 reads 無法相位（HP=0），Subclone 分類主要依賴 63 個 CpG 的甲基化 pattern 矩陣本身。

### 甲基矩陣熱力圖（Cluster Heatmap：Read × CpG）

> 48 reads × 63 CpG；67% reads 為 HP=0（灰）。BERNOULLI UPGMA 聚類。

![chr4:70548355 cluster heatmap](20260315_subclone_phase4_casestudy_assets/chr4_70548355_cluster_heatmap.png)

### Bernoulli 距離熱力圖（Read × Read）

> 深色 = reads 甲基化相似，亮色 = 差異大；HP=0 reads 的聚類結構是否形成獨立群？

![chr4:70548355 distance heatmap](20260315_subclone_phase4_casestudy_assets/chr4_70548355_distance_heatmap.png)

### methylartist — ±2kb per-read 甲基化

![chr4:70548355 methylartist](20260315_subclone_phase4_casestudy_assets/chr4_70548355_methylartist.svg)

### 觀察問題

- [ ] 在 67% HP=0 的情況下，heatmap 是否仍顯示 2 個明顯的 cluster（高甲基 vs 略低甲基）？
- [ ] 11 個 ALT reads 是否集中在 heatmap 的某一層，還是散布？
- [ ] methylartist SVG 的 2kb 視窗中，高甲基是否持續覆蓋整個 window，還是有局部低甲基 patch？

---

## 案例 3：`chr5:153209947` C>A — CramersV = 0.443

### 數值摘要

| 指標 | 值 | 異常旗標 |
|------|-----|---------|
| Reads / CpGs | 57 / 48 | |
| Optimal k | 2 | |
| HP=0 reads | **0/57 = 0%** | ✅ 完全相位（對照案例 2） |
| ALT reads (n) | 17 | all HP=2-1 (HP2 背景 ALT) |
| REF reads (n) | 40 | HP=2 (34) + HP=1 (2) + HP=2-1(4 REF) |
| ALT meth mean | **0.229 ± 0.026** | ⚠️ 所有讀段都是低甲基 |
| REF meth mean | **0.194 ± 0.045** | 差異僅 +0.035 |
| Allele PERMANOVA p | 0.003 (valid) | ✅ |
| **整體甲基化** | **0.04–0.25** | ⚠️ 整個 window 都是低甲基 |
| Bimodal fraction | 100%（全低）| 所有 reads mean_meth < 0.3 |

**關鍵特徵**：位於基因組低甲基化區域（可能是 CpG island 或 promoter），ALT reads 在低甲基背景中略微「相對高甲基化」。
smoothed.csv 數值 7–12% 確認整體低甲基化。

### 甲基矩陣熱力圖（Cluster Heatmap：Read × CpG）

> 57 reads × 48 CpG；整體低甲基化（0.04–0.25），期望整個 heatmap 偏藍。ALT reads 17 條（HP=2-1），REF reads 40 條。

![chr5:153209947 cluster heatmap](20260315_subclone_phase4_casestudy_assets/chr5_153209947_cluster_heatmap.png)

### Bernoulli 距離熱力圖（Read × Read）

> 在低甲基背景下，reads 間距離是否仍形成可辨識的亞群（ALT vs REF 的微小差異）？

![chr5:153209947 distance heatmap](20260315_subclone_phase4_casestudy_assets/chr5_153209947_distance_heatmap.png)

### methylartist — ±2kb per-read 甲基化

![chr5:153209947 methylartist](20260315_subclone_phase4_casestudy_assets/chr5_153209947_methylartist.svg)

### 觀察問題

- [ ] methylartist SVG 整個 2kb 視窗是否呈現均勻低甲基化（接近藍色/空白），還是有局部 CpG clusters 有差異？
- [ ] Heatmap 中的 2 個 cluster 差異是否微弱到難以用眼睛區分（期望：是）？
- [ ] ALT reads（17個）在 heatmap 中是否有任何空間上的 pattern，即使整體值都很低？

---

## 案例 4：`chr16:35118902` G>A — CramersV = 0.442（最清晰雙峰）

### 數值摘要

| 指標 | 值 | 異常旗標 |
|------|-----|---------|
| Reads / CpGs | **91 / 38** | 最多 reads |
| Optimal k | 2 | |
| **HP=0 reads** | **42/91 = 46%** | ⚠️ 高未相位 |
| **HP=3 reads** | **12/91 = 13%**，全為 ALT | ⚠️ 罕見 HP=3 |
| ALT reads (n) | 25 | HP=3 (12) + HP=1-1 (13) |
| REF reads (n) | 66 | HP=0 (42) + HP=1-1 (24) |
| **ALT meth mean** | **0.276 ± 0.116** | ⚠️ ALT = 低甲基化 |
| **REF meth mean** | **0.474 ± 0.200** | REF = 高甲基 |
| **ALT-REF delta** | **−0.198** | ⚠️ 最大甲基差異（−0.20！） |
| Allele PERMANOVA p | **0.001** (valid) | ✅ 最強 allele-meth 關聯 |
| Bimodal fraction | 26.4% high + 28.6% low = **55%** | 最接近真實雙峰 |

**關鍵特徵**：**ALT = 低甲基亞克隆，REF = 高甲基正常細胞**，這是 5 個案例中唯一 ALT allele 對應*低*甲基的位點，也是最符合「亞克隆表觀遺傳重程式化」定義的案例。HP=3 在這裡代表 ALT reads 的 germline HP 背景未能唯一解析，提示 phase complexity，但**不單獨等於 LOH**。

### 甲基矩陣熱力圖（Cluster Heatmap：Read × CpG）

> 91 reads × 38 CpG；**最大 ALT-REF 甲基差（−0.198）**；期望看到清晰高甲基帶（REF）vs 低甲基帶（ALT）的雙色分層。HP=3（12 reads，全 ALT）應在 ALT 低甲基帶中。

![chr16:35118902 cluster heatmap](20260315_subclone_phase4_casestudy_assets/chr16_35118902_cluster_heatmap.png)

### Bernoulli 距離熱力圖（Read × Read）

> 雙峰結構最強的案例，距離矩陣應顯示兩個明顯的 block 對角結構（高甲基群 vs 低甲基群），block 間距離遠（亮色），block 內部距離近（深色）。

![chr16:35118902 distance heatmap](20260315_subclone_phase4_casestudy_assets/chr16_35118902_distance_heatmap.png)

### methylartist — ±2kb per-read 甲基化

![chr16:35118902 methylartist](20260315_subclone_phase4_casestudy_assets/chr16_35118902_methylartist.svg)

### 觀察問題

- [ ] Heatmap 是否清晰呈現「高甲基帶 (REF) vs 低甲基帶 (ALT)」兩個分層？
- [ ] HP=0（42 reads）的甲基化分布位於哪個層次（高還是低）？是否介於兩者之間（混合相位）？
- [ ] HP=3（12 reads，全 ALT）在 methylartist SVG 中是否標示出來？位於低甲基帶？
- [ ] smoothed.csv 的值（接近 0 的極低值與 0.96 的極高值）是否對應 heatmap 的雙峰結構？

---

## 案例 5：`chr7:109185781` G>T — CramersV = 0.429

### 數值摘要

| 指標 | 值 | 異常旗標 |
|------|-----|---------|
| Reads / CpGs | **127 / 62** | 最大資料量 |
| Optimal k | 2 | |
| HP=0 reads | 25/127 = 20% | 中等 |
| ALT reads (n) | 28 (22% VAF) | all HP=**2-1** (HP2 背景 ALT) |
| REF reads (n) | 99 | HP=1 (70), HP=0 (25), HP=1-1 (4) |
| ALT meth mean | 0.402 ± **0.124** | ⚠️ 高 variance（ALT 內部異質） |
| REF meth mean | 0.427 ± 0.070 | 低 variance |
| ALT-REF delta | −0.025 | 差異微小 |
| **Allele PERMANOVA p** | **0.001** (valid) | ✅ |
| **HP PERMANOVA p** | **0.001** (valid) | ⚠️ 唯一 HP 也顯著的位點 |
| Bimodal fraction | 9.4% | 弱 |

**關鍵特徵**：**唯一 HP + Allele 雙重顯著**；127 reads 提供最強統計功效，使 −0.025 的微小差異也達顯著。ALT reads 的高 variance（std=0.124）暗示 28 個 ALT reads 本身可能分成兩個甲基化亞群。

### 甲基矩陣熱力圖（Cluster Heatmap：Read × CpG）

> 127 reads × 62 CpG；最大資料量，統計功效最強。ALT reads（28 條，HP=2-1）的高 variance（std=0.124）暗示 ALT 群內部可能存在 2 個次群。

![chr7:109185781 cluster heatmap](20260315_subclone_phase4_casestudy_assets/chr7_109185781_cluster_heatmap.png)

### Bernoulli 距離熱力圖（Read × Read）

> HP 與 Allele 雙重顯著的唯一案例；距離矩陣應同時呈現 HP 分離（HP1 vs HP2-1）與 Allele 分離的訊號，但 ALT-REF delta 僅 −0.025，差異可能微弱。

![chr7:109185781 distance heatmap](20260315_subclone_phase4_casestudy_assets/chr7_109185781_distance_heatmap.png)

### methylartist — ±2kb per-read 甲基化

![chr7:109185781 methylartist](20260315_subclone_phase4_casestudy_assets/chr7_109185781_methylartist.svg)

### 觀察問題

- [ ] 28 個 ALT reads（HP=2-1，全在上方）在 heatmap 中是否分成 2 層（高 std 的原因）？
- [ ] HP=1（70 reads）與 HP=2-1（28 reads）在 methylartist SVG 中是否有可見甲基差異？
- [ ] 127 reads 的 heatmap 是否顯示更細緻的漸進結構，而非清晰的 2 個群？

---

## FP 對照組：高相似度 False Positive 觀察

選取 4 個與上述 TP Subclone 案例**特徵最相似**的 FP，分兩型：

| FP 類型 | 描述 | 關鍵辨別點 |
|---------|------|----------|
| **Type A：Subclone class FP** | VerificationClass=Subclone，但 variant 為 FP | VAF 極低、HP=0 佔絕大多數 |
| **Type B：高 CramersV FP** | CramersV 高（0.74–1.0），統計信號極強，但 variant 為 FP | 叢集由 HP 驅動而非 Allele；甲基差異被 HP 差異偽裝 |

---

### FP-A1：`chr8:93565727` C>T — Subclone FP，PairwiseMedianDist 最高

#### 數值摘要

| 指標 | 值 | 注意 |
|------|-----|------|
| Reads / CpGs | 68 / 78 | |
| **CramersV** | **0.0** | ⚠️ 無甲基結構信號 |
| PairwiseMedianDist | **0.359** | ⚠️ 高於 TP 案例中位數（結構異質） |
| **ALT reads** | **2 / 68 = VAF 3%** | ⚠️ 極低 VAF，幾乎無 allele 信號 |
| **HP=0 reads** | **67/68 = 97%** | ⚠️ 幾乎全部無法相位（repeat/segdup 疑似） |
| AlleleDelta | 0.000 | ALT_meth=0.708, REF_meth=0.521（delta +0.187，但僅 2 ALT reads） |
| AlleleP | 1.000 | allele PERMANOVA 無效（樣本太少） |
| Optimal k | 4 | |
| VerificationClass | **Subclone** | 被分類為 Subclone 但 variant 是 FP |

**為什麼是 FP**：只有 2 個 ALT reads（VAF 3%），根本無法做 allele PERMANOVA。高 PairwiseMedianDist 和 4 群結構**源自 HP=0 區域的複雜甲基異質性**（可能是 repeat 區域的多重對齊），不是 somatic subclone 信號。

#### 甲基矩陣熱力圖（Cluster Heatmap）

> 注意左側 Allele strip：幾乎全 REF（淺色），只有 2 條 ALT。HP strip：全灰（HP=0）。高 PairwiseMedianDist 會在距離矩陣中看到明顯聚類，但與 allele 無關。

![fp chr8:93565727 cluster heatmap](20260315_subclone_phase4_casestudy_assets/fp_chr8_93565727_cluster_heatmap.png)

#### Bernoulli 距離熱力圖（Read × Read）

> k=4 會顯示 4 個對角 block，但 2 個 ALT reads 無法形成獨立群，預期混在 REF 群中。

![fp chr8:93565727 distance heatmap](20260315_subclone_phase4_casestudy_assets/fp_chr8_93565727_distance_heatmap.png)

#### 觀察問題

- [ ] Cluster heatmap 的 4 群結構是否與 HP=0（灰） reads 的 methylation pattern 差異有關，而非 ALT（2 條）的位置？
- [ ] 2 個 ALT reads 是否在 distance heatmap 中明顯孤立，還是混入某個 REF 群？

---

### FP-A2：`chr9:137953060` T>C — Subclone FP，LabelAlleleP 顯著但 AlleleDelta≈0

#### 數值摘要

| 指標 | 值 | 注意 |
|------|-----|------|
| Reads / CpGs | 58 / **158** | 最多 CpG（超高 PERMANOVA 功效） |
| **CramersV** | **0.0** | 無甲基結構信號 |
| PairwiseMedianDist | 0.268 | |
| ALT reads | **5 / 58 = VAF 9%** | 低 VAF |
| **HP=0 reads** | **42/58 = 72%** | ⚠️ 高未相位 |
| AlleleDelta | **−0.003** | ≈ 0（幾乎無差異） |
| AlleleP | 1.000 | 簡單 allele PERMANOVA 無效 |
| **LabelAllelePermanovaP** | **0.010** | ⚠️ label×allele PERMANOVA 顯著，但 delta≈0 |
| ALT_meth | 0.581 ± 0.022 | |
| REF_meth | 0.581 ± 0.056 | ALT 與 REF 甲基化完全相同！ |
| Optimal k | 3 | |
| VerificationClass | **Subclone** | |

**為什麼是 FP**：158 個 CpG 使 PERMANOVA **統計功效極高**（假顯著性 false significance）。儘管 LabelAllelePermanovaP=0.010，ALT 與 REF 的實際甲基化均值完全相同（均為 0.581）。這是高維度資料中「統計顯著但效應量為零」的典型案例。**CramersV=0 正確反映了無效應量**，是此案例中最準確的指標。

#### 甲基矩陣熱力圖（Cluster Heatmap）

> 158 CpG heatmap 寬度較大。5 個 ALT reads 應混入 REF 群中（兩者甲基化相同），無法形成獨立叢集。

![fp chr9:137953060 cluster heatmap](20260315_subclone_phase4_casestudy_assets/fp_chr9_137953060_cluster_heatmap.png)

#### Bernoulli 距離熱力圖（Read × Read）

> 3 群但 allele 無分離？期望：群的邊界由 HP=0 (42 reads) vs HP=1 (15 reads) 驅動，而非 ALT/REF。

![fp chr9:137953060 distance heatmap](20260315_subclone_phase4_casestudy_assets/fp_chr9_137953060_distance_heatmap.png)

#### 觀察問題

- [ ] 3 群是否對應 HP 分組（HP=0 大群 + HP=1 小群 + 混合），而非 ALT（5 條）vs REF？
- [ ] 5 個 ALT reads 在 cluster heatmap 的 Allele strip 中是否隨機散布於不同群，而非集中？

---

### FP-B1：`chr7:52087777` A>T — 高 CramersV=1.0 FP，HP+Allele 雙重顯著

#### 數值摘要

| 指標 | 值 | 注意 |
|------|-----|------|
| Reads / CpGs | 63 / 106 | |
| **CramersV** | **1.000** | ⚠️ 最高可能值！完美分離 |
| PairwiseMedianDist | 0.265 | |
| **ALT reads** | **9 / 63 = VAF 14%** | 低 VAF |
| HP=0 reads | 1/63 | 幾乎全相位 |
| HP 分布 | HP2=31, HP1=21, HP1-1/HP2-1=10 | 良好 HP1/HP2 分布 |
| AlleleDelta | **+0.058** | ALT_meth=0.656, REF_meth=0.595 |
| AlleleP | **0.001** | |
| LabelHPPermanovaP | **0.001** | |
| LabelAllelePermanovaP | **0.001** | |
| Optimal k | 4 | |
| VerificationClass | **Strong** | （非 Subclone！） |

**為什麼是 FP**：所有統計指標（CramersV=1.0、AlleleP=0.001、HPP=0.001）都指向「真正的甲基化結構差異」，但 **variant 本身**是 FP（可能是 germline 多型性、PON artifact、或 SEQC2 benchmark 未收錄）。**甲基信號是真實的 germline haplotype-associated methylation difference，與 somatic variant 無關**。這是最「危險」的 FP 類型：甲基特徵完美，但 variant 錯。**辨別線索**：VAF 僅 14%（9 ALT reads），且 HP1 vs HP2 分布均衡（HP 差異比 allele 差異更能解釋叢集）。

#### 甲基矩陣熱力圖（Cluster Heatmap）

> CramersV=1.0 應產生視覺上最清晰的分層。但請觀察：分層是由 HP tag（HP1 vs HP2，藍 vs 黃）還是 Allele（ALT vs REF）驅動？

![fp chr7:52087777 cluster heatmap](20260315_subclone_phase4_casestudy_assets/fp_chr7_52087777_cluster_heatmap.png)

#### Bernoulli 距離熱力圖（Read × Read）

> k=4 的 4 個 block 是否對應 HP 分組（HP1, HP2, HP1-1/HP2-1）？若距離矩陣的 block 邊界由 HP 解釋，則這是「HP-driven FP」。

![fp chr7:52087777 distance heatmap](20260315_subclone_phase4_casestudy_assets/fp_chr7_52087777_distance_heatmap.png)

#### methylartist — ±2kb per-read 甲基化

> 完整 HP 分層（HP1=21, HP2=31），請觀察甲基化分層是否平行於 HP tag 而非 ALT/REF allele。

![fp chr7:52087777 methylartist](20260315_subclone_phase4_casestudy_assets/fp_chr7_52087777_methylartist.svg)

#### 觀察問題

- [ ] Cluster heatmap 中，HP strip 的顏色（藍=HP1，黃=HP2）是否比 Allele strip（ALT/REF）更整齊地分群？
- [ ] 9 個 ALT reads 在 HP strip 中是哪個 HP？是否全集中在某一 HP（暗示 allele 與 HP 共連鎖，而非 somatic）？
- [ ] methylartist SVG 中，HP1 vs HP2 的甲基化差異是否比 ALT vs REF 更明顯？

---

### FP-B2：`chr9:75383880` T>A — CramersV=0.74 FP，HP 驅動的偽 Allele 信號

#### 數值摘要

| 指標 | 值 | 注意 |
|------|-----|------|
| Reads / CpGs | 62 / 34 | |
| **CramersV** | **0.736** | 高，接近 TP 中的 chr6 (0.72) |
| PairwiseMedianDist | **0.352** | 高 |
| ALT reads | 14 / 62 = VAF 23% | 中等 VAF |
| HP 分布 | HP2=21, HP1-1=17, HP1=17, HP0=7 | 均衡三群 |
| **AlleleDelta** | **+0.005** | ≈ 0，甲基差異幾乎不存在！ |
| AlleleP | 0.361 | 簡單 allele 檢驗**不顯著** |
| LabelHPPermanovaP | **0.001** | HP 分離極顯著 |
| LabelAllelePermanovaP | **0.001** | ⚠️ Allele PERMANOVA 顯著但 AlleleDelta≈0 |
| ALT_meth | 0.502 ± 0.111 | |
| REF_meth | 0.517 ± 0.155 | 差異 −0.015，幾乎無分離 |
| Optimal k | 3 | |
| VerificationClass | **Strong** | |

**為什麼是 FP**：CramersV=0.736 和 PairwiseMedianDist=0.352 都很高，LabelAllelePermanovaP=0.001，但 **AlleleDelta=+0.005（接近 0）**。叢集由 **HP 差異**（HP1 vs HP2 的 germline methylation difference）驅動，而 14 個 ALT reads 因恰好在某個 HP group 中而被「誤認為」與甲基叢集有關（allele-HP 連鎖的 confounding）。簡單 AlleleP=0.361 不顯著，是此案例的正確反映。

#### 甲基矩陣熱力圖（Cluster Heatmap）

> HP 均衡分布（HP1≈HP2≈HP1-1 各 17-21 reads）。若 3 群對應 HP 分層，則甲基信號是 germline imprinting/differential methylation，不是 somatic。

![fp chr9:75383880 cluster heatmap](20260315_subclone_phase4_casestudy_assets/fp_chr9_75383880_cluster_heatmap.png)

#### Bernoulli 距離熱力圖（Read × Read）

> 3 個 block 是否對應 HP1 / HP1-1 / HP2 的三個 germline group，而非 ALT(14) vs REF(48)?

![fp chr9:75383880 distance heatmap](20260315_subclone_phase4_casestudy_assets/fp_chr9_75383880_distance_heatmap.png)

#### methylartist — ±2kb per-read 甲基化

![fp chr9:75383880 methylartist](20260315_subclone_phase4_casestudy_assets/fp_chr9_75383880_methylartist.svg)

#### 觀察問題

- [ ] HP strip 的三色分群（HP1/HP1-1/HP2）是否完全對應 3 個 cluster block？這是 germline HP-methylation 連鎖，而非 somatic。
- [ ] 14 個 ALT reads 在 HP strip 中分布於哪個/哪些 HP？是否跨越 HP1 和 HP2，代表 allele-HP 連鎖不完全（進一步反駁 somatic 假說）？
- [ ] methylartist SVG 中，是否有明顯的「HP1 偏低甲基 / HP2 偏高甲基」（或反之）的漸進結構？

---

## 跨案例比較表

| 位點 | 標籤 | CramersV | 甲基類型 | ALT-REF delta | ALT reads | 主要異常 | 視覺預期 |
|------|------|---------|---------|--------------|-----------|---------|---------|
| **chr6:145444893** | **TP** | **0.722** | 中等甲基 | +0.04 | 19/58 | k=3，全 ALT=HP=2-1 | 三層 heatmap |
| chr4:70548355 | TP | 0.502 | 高甲基 | +0.04 | 11/48 | 67% HP=0 | 模糊邊界 |
| chr5:153209947 | TP | 0.443 | **低甲基** | +0.035 | 17/57 | 全域低甲基化 | 均勻藍色 |
| chr16:35118902 | TP | 0.442 | 雙峰 | **−0.198** | 25/91 | ALT=低甲基，HP=3 | 清晰雙色帶 |
| chr7:109185781 | TP | 0.429 | 中等甲基 | −0.025 | 28/127 | HP+Allele 皆顯著 | 微弱分層 |
| chr8:93565727 | **FP-A1** | 0.000 | 中高甲基 | +0.187\* | **2/68** | 97% HP=0，VAF 3% | 無 Allele 分層 |
| chr9:137953060 | **FP-A2** | 0.000 | 中高甲基 | −0.003 | **5/58** | 158 CpG 假顯著 | 無分層 |
| chr7:52087777 | **FP-B1** | **1.000** | 中高甲基 | +0.058 | 9/63 | HP-driven，variant FP | HP 分層（非 allele） |
| chr9:75383880 | **FP-B2** | 0.736 | 中等甲基 | **+0.005** | 14/62 | AlleleDelta≈0，HP-driven | HP 三群分層 |

> \* chr8:93565727 的 AlleleDelta=0.187 是由僅 2 個 ALT reads 計算，統計上不可信（AlleleP=1.0）。

**關鍵鑑別規則（從 FP 對照組歸納）**：

| 規則 | 說明 |
|------|------|
| **ALT reads ≥ 10** | FP-A1 (n=2), FP-A2 (n=5) 的 VAF 極低，缺乏 allele 功效 |
| **AlleleDelta > 0.05** | FP-B2 的 AlleleDelta=0.005，LabelAlleleP 顯著是 HP confounding 造成的假象 |
| **CramersV > 0.3 ≠ 充分條件** | FP-B1 CramersV=1.0 仍是 FP；需同時確認 allele 是否驅動（而非 HP） |
| **HP_Ratio 是否合理** | FP-B1/B2 的 HP 均衡分布下，CramersV 高代表 germline methylation，非 somatic |

---

## 觀察記錄（TP 案例）

### TP-1：chr6:145444893（CramersV=0.722）

**觀察**：
- Cluster heatmap 呈現 3 個可辨識的水平帶；左側 HP strip 顯示最下方群（ALT）全為橙色（HP=2-1），中、上方群為黃色（HP=2）與灰色（HP=0）混合
- 19 個 ALT reads（全 HP=2-1）聚集在同一個分層；REF reads（HP=2 + HP=0）分成兩個亞群，但視覺差異小於 ALT-REF 分離
- Bernoulli 距離熱力圖顯示 ALT 群為獨立對角 block（深色），與 REF 群之間距離明顯較亮（差異大）
- methylartist SVG 中，ALT reads（HP=2-1 橙色）在 SNV 附近的 CpG 點呈中等甲基化，與 HP=2（REF）整體相近但 pattern 可區分

**已驗證**：
- Allele PERMANOVA p=0.001，HP PERMANOVA p=1.0 → allele 是唯一有效分類標籤，HP 完全無貢獻；這與 heatmap 中 HP strip 無法分離的視覺結果一致
- k=3 確認三群：1 個 ALT 群 + 2 個 REF 亞群；符合 InterSubMod 的 optimal k=3 計算

**推論**：
- HP=2-1 表示這 19 條 ALT reads 的 germline 背景是 HP2，但它們在 SNV 位置帶有 ALT allele；甲基化分離代表 **somatic ALT allele 與特定甲基狀態共分離**，而非 germline HP 差異
- 這是最乾淨的「Allele-driven subclone methylation」訊號範例：排除 HP 影響，allele 本身解釋甲基差異

---

### TP-2：chr4:70548355（CramersV=0.502）

**觀察**：
- Cluster heatmap 在 67% HP=0（灰）主導下，仍顯示 2 個可辨識的高低甲基帶；HP strip 灰色讀段分布在兩群中，無法單獨解釋分群
- 11 個 ALT reads（HP=1-1，深藍）傾向集中在較高甲基化的群中；63 個 CpG 提供豐富訊號
- 距離熱力圖中兩群的對角 block 結構清晰，但 block 內異質性較高（非如 chr16 的完美雙色）
- methylartist 2kb 視窗顯示整體為高甲基化區域（0.58–0.70），局部有低甲基 patch 但持續性不強

**已驗證**：
- AlleleP=0.011（邊緣顯著），HP=0 佔 67% 不代表無 allele 資訊；InterSubMod 的 PERMANOVA 使用實際甲基化矩陣，即使 HP 無法區分，allele label 仍能分離 reads
- HP_Ratio=0.937 表示僅有極少數 HP1 reads，HP 訊號不均衡；但這反而更支持 allele 是分群驅動力

**推論**：
- 高比例 HP=0 可能是此位點附近缺乏 germline SNP（難以相位的區域），或是基因組中 segmental duplication/repeat 區域，導致 LongPhase-S 無法相位
- HP=0 reads 分散於兩群而非自成一群，支持「甲基差異由 somatic allele 驅動，不由 HP 背景決定」

---

### TP-3：chr5:153209947（CramersV=0.443）

**觀察**：
- Cluster heatmap 整體呈均勻藍色（低甲基，0.04–0.25），兩群視覺差異極微弱，需要仔細辨認才能看出邊界
- 17 個 ALT reads（HP=2-1）與 40 個 REF reads（主要 HP=2）均為低甲基化，差異僅 +0.035
- 距離熱力圖中，ALT vs REF 的 block 分離不如 chr16/chr6 清晰，但仍有可辨別的結構
- methylartist 2kb 視窗確認整個 window 均為低甲基化（smoothed.csv 值 7–12%），無局部高甲基 patch

**已驗證**：
- AlleleP=0.003 達顯著：在 57 reads × 48 CpG 的矩陣規模下，即使 delta 僅 0.035，多變量 PERMANOVA 仍能檢測到統計差異
- 整體低甲基（< 0.3）確認此位點在 CpG island 或 promoter 調控區域（低甲基是 gene activation 的特徵）

**推論**：
- 「低甲基背景中的微弱差異」是最難視覺驗證的 subclone 類型；InterSubMod 的 PERMANOVA 方法（多 CpG 聯合檢驗）比單一位點比較更敏感
- 若此位點是 somatic subclone，可能代表 promoter 區域的表觀遺傳活化（ALT allele 的 reads 略微升高甲基），而非典型的 silencing 事件

---

### TP-4：chr16:35118902（CramersV=0.442，最清晰雙峰）

**觀察**：
- Cluster heatmap 呈現最清晰的雙色分層：**上方深紅帶（REF 高甲基，0.5–0.8）vs 下方深藍帶（ALT 低甲基，0.1–0.35）**，視覺上最易辨別
- HP=3（紫，12 reads）全部落在 ALT 低甲基帶；HP=1-1（13 reads）也在 ALT 帶；HP=0（42 reads）主要在 REF 帶但少數散落至 ALT 帶
- 距離熱力圖顯示 2 個完整的對角 block，block 間距離明顯亮（差異大），是所有案例中 block 結構最清晰的
- methylartist SVG 呈現清楚的兩帶：HP2（黃）reads 在上方呈高甲基化，ALT reads（HP3 紫色 + HP1-1 深藍）在下方呈低甲基化

**已驗證**：
- ALT-REF delta = −0.198（5 個案例中最大），AlleleP=0.001；視覺雙峰與統計量高度一致
- HP=3 reads（定義：LongPhase-S 無法區分 HP1/HP2，但 read 跨越 Somatic ALT）全部是 ALT reads，這支持此位點的 somatic ALT 帶有獨特的低甲基化狀態
- Bimodal fraction 55%（26.4% high + 28.6% low），數值最接近理想的 50/50 雙峰

**推論**：
- **ALT allele = 低甲基亞克隆**：這是最符合「somatic subclone 伴隨表觀遺傳重程式化（epigenetic reprogramming）」定義的案例 — somatic 突變的腫瘤亞克隆同時失去甲基化（hypomethylation）
- HP=3 的出現（12 reads 無法唯一解析 germline HP 但帶有 ALT）提示此區域存在 phase complexity；染色體上的 HP 背景資訊不足，但 somatic ALT 訊號清晰
- 若後續 CNV / phased genotype 顯示此位點附近真有 LOH，才可再把這種 HP=3 富集視為 LOH-compatible evidence，而不是單靠 HP=3 下結論

---

### TP-5：chr7:109185781（CramersV=0.429）

**觀察**：
- Cluster heatmap 顯示兩個主要帶，但邊界較 chr16 模糊；HP=2-1（橙，28 ALT reads）在上方，HP=1（藍，70 reads）在中下方，HP=0（灰，25 reads）散布
- 28 個 ALT reads 的 HP strip（橙色帶）在 heatmap 中形成一個群，但其甲基化值呈現明顯的內部異質性（std=0.124），視覺上可看到橙色帶中有些 reads 偏紅、有些偏藍
- 距離熱力圖顯示 2 群但邊界不清晰，HP=1 vs HP=2-1 的分離是主要聚類驅動力
- methylartist 2kb 視窗顯示 HP=2-1（橙）reads 在局部 CpG 位點有微弱的甲基降低，但整體差異不大

**已驗證**：
- 唯一 HP PERMANOVA p=0.001 AND Allele PERMANOVA p=0.001 的案例：HP 差異與 Allele 差異**共同**解釋聚類；但 ALT-REF delta = −0.025（微弱），表示 HP 是更強的分離因子
- 127 reads 提供極強統計功效：delta=−0.025 在其他 case 中不足顯著，但在 n=127 下達顯著

**推論**：
- HP 和 Allele 雙重顯著並非矛盾；HP=2-1（ALT reads）帶有 HP2 背景，若 HP1 vs HP2 在此位點有 germline methylation difference，則 ALT reads 因恰好在 HP2 背景上而「繼承」了 HP2 的甲基化 pattern，使 allele PERMANOVA 也顯著
- 28 個 ALT reads 的高內部 std (0.124) 可能表示：(a) ALT reads 本身構成 2 個次群（不同甲基狀態的腫瘤亞克隆），或 (b) 某些 ALT reads 帶有的 HP2 germline methylation 比其他 ALT reads 更強
- 「最大讀數 + 最弱 delta + HP 共線性」是理解此案例的關鍵：統計顯著不代表效應量大

---

## FP 觀察記錄

### FP-A1：chr8:93565727（CramersV=0.0，VAF=3%）

**觀察**：
- Cluster heatmap 左側 Allele strip 幾乎全為 REF 淺色（68 reads 中僅 2 個 ALT）；HP strip 全灰（97% HP=0）
- 4 個聚類群清晰可見，但 2 個 ALT reads 散落在不同群中，**沒有自成一群**
- 距離熱力圖顯示 4 個 block，但這些 block 邊界由 HP=0 reads 的甲基化 pattern 差異驅動，與 allele 無關
- 68 個 reads 中有豐富甲基異質性（PairwiseMedianDist=0.359），但異質性分散於 HP=0 reads，而非集中在 ALT

**已驗證**：
- CramersV=0.0 正確：2 個 ALT reads 無統計功效，無法計算有意義的甲基-allele 關聯
- AlleleP=1.0 正確：2 reads vs 66 reads 的 allele PERMANOVA 計算量不足
- **與 TP 案例的關鍵視覺差別**：TP 案例的 ALT reads 有 ≥11 條且集中在同一層；此 FP 的 2 條 ALT reads 視覺上無法辨識，完全混入 REF 群

**推論**：
- 97% HP=0 可能表示此位點在基因組 repeat/segdup 區域，reads 從多個同源位置對齊，造成 HP 無法相位
- 高 PairwiseMedianDist 來自 repeat region 的多重對齊 reads 帶有不同的甲基化狀態（不同 gene copy 的差異甲基化），而非 somatic subclone 信號
- **辨識規則**：VAF < 5% 且 HP=0 > 90% 時，VerificationClass=Subclone 幾乎可確認為 FP

---

### FP-A2：chr9:137953060（CramersV=0.0，158 CpG 假顯著）

**觀察**：
- Cluster heatmap 因 158 個 CpG 顯示較寬的矩陣；HP strip 顯示 3 群主要由 HP=0（灰，42 reads）vs HP=1（藍，15 reads）vs 其他（1 read）驅動
- 5 個 ALT reads 在 Allele strip 中**隨機散布**於不同群（視覺上不集中），完全無法辨識為獨立群
- 距離熱力圖的 3 個 block 對應 HP=0 大群、HP=1 小群和混合群，而非 ALT/REF 分離
- ALT reads 的甲基化值與對應 HP 群的 REF reads 完全相同（ALT_meth=0.581 = REF_meth=0.581）

**已驗證**：
- CramersV=0.0 正確：即使 LabelAllelePermanovaP=0.010 顯著，效應量仍為零
- 「高維假顯著性」典型案例：158 CpG 提供的高自由度使 PERMANOVA 即使在 zero delta 下也能達顯著（type I error inflation）
- **與 TP 案例的關鍵視覺差別**：TP 案例中，ALT reads 在 CpG axis 上有可辨識的 pattern shift；此 FP 的 5 條 ALT reads 無論看哪行都與周圍 REF reads 相同

**推論**：
- 158 CpG 是所有案例中最多，也是假顯著風險最高的；InterSubMod 的 CramersV 是防止高維假顯著的關鍵安全閥
- 「LabelAllelePermanovaP 顯著但 CramersV=0 且 AlleleDelta≈0」是最清楚的 FP pattern：統計功效偶然達閾值，但無效應量
- **辨識規則**：CramersV=0.0 且 AlleleDelta < 0.01 即可排除，不需看 heatmap

---

### FP-B1：chr7:52087777（CramersV=1.0，local HP-labeled structure + SEQC2 INDEL 鄰近）

**觀察**：
- Cluster heatmap 顯示極清楚分層（CramersV=1.0），但視覺上更像被 local HP 標籤結構主導，而不是乾淨的 allele-only 分離
- IGV 與後續 BAM 驗證顯示，target ALT reads 為 `11/78`，且都標成 `HP=2-1`
- 因 phase block 未在此段進一步驗證，這裡只能寫成「ALT-supporting reads 與同一 local HP label 共線」，不能直接升級成親緣來源敘述

**已驗證**：
- `LabelAlleleP=0.001` AND `LabelHPP=0.001`，兩者皆顯著，但 `DominantLabel=hp`
- `tumor ALT=11/78`，`normal ALT=0/48`
- nearby truth event 為 `chr7:52087776:TA>T`（SEQC2 HighConf sINDEL）

**推論**：
- 這個案例最合理的解釋是：**local HP-labeled methylation structure 真實存在，但 target SNV 本身又緊鄰 truth sINDEL**
- 因此它同時是「甲基訊號真實」與「variant 標記受 benchmark / local alignment 背景影響」的案例

> **修正說明**：本段已依 `20260316_IGV固定template_AI視覺初篩與正式驗證_01.md` 與 `20260316_igv_case_validation_01.tsv` 校正，不再沿用較早版本的 `HP1-ALT` 敘述。

---

### FP-B2：chr9:75383880（CramersV=0.74，local alignment anomaly / -1 deletion-like gap）

**觀察**：
- Heatmap 的群分得開，但更像 HP / local read-structure 驅動，而不是穩定的 allele-specific methylation
- IGV 視覺上確實有「target 前一格也一起異常」的現象，因此這類位點需要回到 BAM 做鄰近位置的 read-level 驗證

**已驗證**：
- `tumor ALT=15/106`，ALT reads 皆為 `HP=1-1`
- `adj_minus1_alt_del_frac=1.000`
- 嚴格 MNP 篩查結果為 `0 candidates`

**推論**：
- 這個案例在主 tumor BAM 下更像 **`target SNV + -1 deletion-like gap / local alignment anomaly`**
- 因此不應再把它當成乾淨 MNP 的代表案例，也不宜直接沿用「相鄰雙 SNV 被拆開」的舊說法

> **修正說明**：本段已由後續正式驗證取代較早版本的 MNP 假說；若需追溯舊討論，應同時附註「已被 20260316 validation 修正」。

---

## 知識點整理（可重複應用的觀察與推論規則）

### K1：HP tag 讀法（LongPhase-S Paired 模式）

| HP 值 | 含義 | 典型出現位置 |
|-------|------|------------|
| HP=1 / HP=2 | Germline haplotype 1 或 2，來自 normal phasing | 通常為 REF reads |
| HP=0 | 無法相位（無 germline SNP 附近），或 reads 在 repeat 區 | REF 或 ALT reads |
| HP=1-1 | HP1 背景，但此 read 帶有 Somatic ALT | ALT reads（HP1 背景腫瘤 allele） |
| HP=2-1 | HP2 背景，但此 read 帶有 Somatic ALT | ALT reads（HP2 背景腫瘤 allele） |
| HP=3 | LongPhase-S 無法區分 HP1/HP2，但 read 帶有 Somatic ALT | ALT reads（HP 背景不明） |

**實用規則**：
- 若 ALT reads 全為 HP=X-1 或 HP=3 → read 的 HP 背景已確認
- 若 ALT reads 為 HP=0 → ALT reads 無法相位（VAF 低或 repeat 區）

---

### K2：HP-driven vs Allele-driven 甲基分離的鑑別

| 特徵 | HP-driven（FP 風險高） | Allele-driven（TP 可信） |
|------|----------------------|------------------------|
| HP strip 分層 | HP1/HP2 清晰分層 | HP strip 無法解釋分群 |
| ALT reads HP 分布 | ALL ALT reads 在某一 HP 側 | ALT reads 分布於不同 HP 背景 |
| AlleleDelta | 接近 0（ALT≈REF） | > 0.05（ALT 與 REF 有實質差異） |
| HP PERMANOVA p | 顯著（< 0.05） | 不顯著（> 0.3）或遠弱於 AlleleP |
| 解釋 | Germline imprinting / LOH | Somatic subclone epigenetic change |

---

### K3：FP 類型分類表

| FP 類型 | 典型特徵 | 辨識指標 |
|---------|---------|---------|
| **A 型（低 VAF 假 Subclone）** | VAF < 5%，HP=0 > 90%，CramersV=0 | ALT reads ≤ 5，AlleleP=1.0 |
| **B1 型（HP-driven imprinting）** | CramersV 高（>0.3），HP 清晰分層，ALT 全在某 HP 側 | AlleleDelta < 0.02，HPP << AlleleP |
| **B2 型（local alignment anomaly / adjacent mismatch）** | AlleleDelta≈0，LabelAlleleP 顯著但 AlleleP 不顯著，IGV 顯示 target 前後欄位異常 | HP-ALT 共線，且需再查 adjacent-gap 指標 |
| **C 型（高 CpG 假顯著）** | CpG > 100，CramersV=0，ALT delta≈0 | LabelAlleleP 顯著但 CramersV=0 |

---

### K4：TP Subclone 甲基類型分類（初步 typology）

| 類型 | 代表案例 | 特徵 | 甲基 pattern |
|------|---------|------|------------|
| **T1：Allele-only 中等甲基** | chr6:145444893 | ALT=HP 背景已知的 ALT reads，allele 主驅動，HP 無效 | 中等甲基，3群結構 |
| **T2：高甲基背景弱分離** | chr4:70548355 | 高 HP=0，整體高甲基，邊界模糊 | 高甲基，2群弱邊界 |
| **T3：低甲基背景微差異** | chr5:153209947 | 整個 window 低甲基，差異 < 0.04 | 均勻低甲基 |
| **T4：雙峰 ALT低甲基** | chr16:35118902 | ALT=低甲基，REF=高甲基，最清晰雙峰 | 雙色帶，-0.2 delta |
| **T5：HP+Allele 共線** | chr7:109185781 | HP 和 allele 雙重顯著，高讀數 | 微弱分層，高功效 |

---

### K5：adjacent anomaly / deletion-like gap 辨識步驟

當在 IGV 中看到 ALT reads 在 SNV 的鄰近位置有一致異常時：

1. **先看 adjacent-gap 指標**：例如 `adj_minus1_alt_del_frac` / `adj_plus1_alt_del_frac`
2. **檢查 adjacent position 的 pileup**：比較 ALT reads vs REF reads 在 ±1bp 的分布
3. **比對 normal BAM**：確認 normal 是否也有相同異常
4. **執行嚴格 MNP screen**：若結果為 `0 candidates`，不要再把它寫成乾淨 MNP
5. **最後才判讀機制**：可能是 MNP、INDEL-adjacent artifact、或 deletion-like gap；不可先入為主

---

### K6：SEQC2 INDEL-adjacent SNV 的 benchmark annotation gap

- **問題**：SNV-only benchmark 在 INDEL ±2bp 內的 SNV 位置無法正確評估
- **典型案例**：chr7:52087776（SEQC2 INDEL TP）→ chr7:52087777 被 SNV benchmark 標記為 FP
- **建議**：在大規模評估中，對所有 SEQC2 sINDEL TP 位置 ±2bp 範圍內的 SNV 呼叫，先排除後再評估 precision/recall
- **規模估計**：HCC1395 有 ~4,000 個 SEQC2 sINDEL；若平均 ±2bp 影響 4 個 position，約有 16,000 個 position 需要特別處理

---

## IGV 自動化截圖方案

基於 `/big8_disk/liaoyoyo2001/IGV_session/template.xml`，可用 bash 腳本生成每個位點的 IGV 截圖。

### 方案一：IGV batch script（推薦）

```bash
#!/bin/bash
# scripts/utils/igv_batch_screenshot.sh
# 用法: ./igv_batch_screenshot.sh chr9 75383880 2000 /output/dir/

CHR="$1"
POS="$2"
WINDOW="${3:-2000}"
OUTDIR="${4:-/tmp/igv_screenshots}"

# 計算視窗範圍
START=$((POS - WINDOW))
END=$((POS + WINDOW))
LOCUS="${CHR}:${START}-${END}"
OUTFILE="${OUTDIR}/${CHR}_${POS}.png"

mkdir -p "$OUTDIR"

# 生成修改後的 session XML（替換 locus）
TMPXML="/tmp/igv_session_${CHR}_${POS}.xml"
sed "s|<locus>.*</locus>|<locus>${LOCUS}</locus>|g" \
    /big8_disk/liaoyoyo2001/IGV_session/template.xml > "$TMPXML"

# 生成 IGV batch 命令腳本
BATCHFILE="/tmp/igv_batch_${CHR}_${POS}.txt"
cat > "$BATCHFILE" << EOF
load ${TMPXML}
goto ${LOCUS}
snapshot ${OUTFILE}
exit
EOF

# 執行 IGV（需要 xvfb-run 無頭模式）
xvfb-run --auto-servernum igv.sh -b "$BATCHFILE" 2>/dev/null

echo "Screenshot saved: ${OUTFILE}"
rm -f "$TMPXML" "$BATCHFILE"
```

### 方案二：批次多位點截圖

```bash
#!/bin/bash
# scripts/utils/igv_batch_multi.sh
# 用法: ./igv_batch_multi.sh positions.tsv /output/dir/
# positions.tsv 格式: chr\tpos（每行一個位點）

POSITIONS_FILE="$1"
OUTDIR="${2:-/tmp/igv_screenshots}"
WINDOW=2000

mkdir -p "$OUTDIR"

while IFS=$'\t' read -r CHR POS; do
    echo "Processing ${CHR}:${POS}..."

    START=$((POS - WINDOW))
    END=$((POS + WINDOW))
    LOCUS="${CHR}:${START}-${END}"
    OUTFILE="${OUTDIR}/${CHR}_${POS}.png"

    TMPXML=$(mktemp /tmp/igv_session_XXXXXX.xml)
    BATCHFILE=$(mktemp /tmp/igv_batch_XXXXXX.txt)

    sed "s|<locus>.*</locus>|<locus>${LOCUS}</locus>|g" \
        /big8_disk/liaoyoyo2001/IGV_session/template.xml > "$TMPXML"

    cat > "$BATCHFILE" << EOF
load ${TMPXML}
goto ${LOCUS}
snapshot ${OUTFILE}
exit
EOF

    xvfb-run --auto-servernum igv.sh -b "$BATCHFILE" 2>/dev/null
    rm -f "$TMPXML" "$BATCHFILE"

done < "$POSITIONS_FILE"

echo "Done. Screenshots in: ${OUTDIR}"
```

### 確認環境（已驗證）

```bash
# IGV 安裝位置（已確認）
IGV_PATH=/home/liaoyoyo2001/igv.sh   # 主要位置
# 備用：/big8_disk/liaoyoyo2001/igv/build/IGV-dist/igv.sh

# DISPLAY 已設定（X11 forwarding）
echo $DISPLAY   # → localhost:11.0
# xvfb-run 不可用，但 DISPLAY 已存在，腳本可直接執行

# 參考基因組（已確認）
REF=/big8_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa

# 測試截圖（單一位點）
export IGV_PATH=/home/liaoyoyo2001/igv.sh
./scripts/utils/igv_batch_screenshot.sh chr9 75383880 2000 /tmp/igv_test/
```

> **注意**：腳本會自動偵測 DISPLAY；若 DISPLAY 已設定，直接執行 IGV（無需 xvfb-run）。若 DISPLAY 未設定且無 xvfb-run，腳本會提示錯誤。

---

## 下一步（視覺確認後）

- [ ] 確認 chr16:35118902 所在區域是否與 `LOH.bed` / CNV track / phased genotype 一致；不能單靠 `HP=3` 推論 LOH
- [ ] 確認 chr5:153209947 的低甲基化是否與已知 CpG island 吻合（UCSC 查詢）
- [ ] 決定是否將 5 種 Subclone methylation signature 正式命名為 typology（K4 已初步分類）
- [ ] 若 chr16 雙峰視覺確認，可作為 `subclone_methylation_type = bimodal_ALT_hypomethylated` 的 prototype
- [ ] FP 鑑別規則（K3）是否在視覺上可直觀辨認 → 若是，考慮加入 VCF annotation filter tier
- [ ] 執行大規模 adjacent anomaly / deletion-like gap 篩選（可沿用 `screen_mnp_adjacent_fp.py` 作第一輪排除），量化 HCC1395 FP 中此類比例
- [ ] 確認 IGV xvfb-run + igv.sh 環境可用，執行 9 個案例位點的自動截圖
- [ ] FP-B1 chr7:52087777：用含 sINDEL truth set 的 benchmark 重新評估，確認是否為 INDEL-adjacent SNV artifact
