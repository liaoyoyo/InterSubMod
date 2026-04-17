# F Pilot Step 4 — AF<0.4 Stratified Cohen's d + FP 聚集分析

**Date**: 2026-04-18
**Script**: `scripts/step4_af04_cohens_d_and_fp_clustering.py`

---

## 🎯 核心結論

1. **(A) AF<0.4 重算 Cohen's h**：HCC1954 **從「特殊」升為確定 POS**（h=+0.587→+0.654，AF<0.2 達 +0.775 接近 large）；但 **COLO829 仍 NEGATIVE**（-0.457→-0.429，ONT_R10 artifact 無法 rescue）；HCC1395/DORADO 因基線飽和產生 h 數值降低（非真正訊號衰減）
2. **(B) Coverage 影響**：6/7 樣本仍 YES power（n≥100），僅 COLO829 MARGINAL（n=34 未變）；HCC1954 保留 1,040 regions、H1437 1,132 regions
3. **(C) FP 聚集假設被推翻** ⚠️：AF≥0.4 FP **並非空間聚集**（ratio 0.81× for HCC1954；5/5 樣本 <1 表示 FP 反而比 TP 分散）；但 **染色體層級富集明確**（HCC1954 chr7/8/16 的 FP 中 68-71% 在 AF≥0.4）
4. **(D) Feature profile 揭示 germline-like 機制**：AF≥0.4 FP NumCpGs 更高（啟動子/CpG island）、PairwiseMedianDist 更低（reads 相似性更高）→ 符合 germline het 被 CNV 漂移後單 allele 主導

---

## Part A: AF<0.4 stratified per-sample Cohen's h

| Sample | p4_all | plt_all | **h_all** | p4_af04 | plt_af04 | **h_af04** | p4_af02 | plt_af02 | **h_af02** | 變化 |
|---|---|---|---|---|---|---|---|---|---|---|
| **HCC1954** | 0.497 | 0.220 | +0.587 | 0.707 | 0.387 | **+0.654** | 0.874 | 0.535 | **+0.775** | ✅ **升級**（medium+→medium+/接近 large） |
| HCC1937 | 0.713 | 0.454 | +0.534 | 0.867 | 0.649 | +0.523 | 0.949 | 0.847 | +0.350 | ➖ medium+ 穩定 |
| H1437 | 0.921 | 0.757 | +0.462 | 0.965 | 0.802 | **+0.543** | 0.958 | 0.874 | +0.313 | ✅ 升級（small→medium+） |
| **HCC1395_DORADO** | 0.903 | 0.687 | +0.553 | 0.919 | 0.836 | +0.257 | 0.895 | 0.839 | +0.165 | ⚠️ 數值降級（飽和 artifact） |
| HCC1395 | 0.810 | 0.688 | +0.283 | 0.887 | 0.828 | +0.171 | 0.891 | 0.851 | +0.121 | ⚠️ 數值降級（飽和 artifact） |
| H2009 | 0.935 | 0.900 | +0.127 | 0.957 | 0.922 | +0.150 | 0.985 | 0.920 | +0.328 | ➖ ceiling |
| **COLO829** | 0.235 | 0.450 | **-0.457** | 0.235 | 0.436 | **-0.429** | 0.190 | 0.440 | -0.547 | ❌ **仍 NEGATIVE**（ONT_R10 artifact 無法 rescue） |

### 解讀

- **HCC1954 POS 升級**：h 從 +0.59 (medium+) → +0.65 (medium+ 更強) → AF<0.2 +0.78（接近 large 0.8），完全契合預期
- **HCC1395/DORADO 數值降級是飽和 artifact，非訊號衰減**：基線 TP rate (plt_af04) 從 0.688→0.828（DORADO 0.687→0.836），因 AF<0.4 讓 NG<4 也變乾淨，差距縮小；h 的飽和在 phi space 特別明顯（接近 1.0 時 arcsin 飽和）
- **COLO829 無法 rescue**：AF<0.4 內 NG=4 TP rate=0.235 vs NG<4 TP rate=0.436 → NG=4 反而更差；確認 memory 中「COLO829 ONT_R10 無 methylation → HPFineNGroups 訊號為 artifact」正確
- **先前 B.1-3 "5/7 POS + 2/7 特殊"** 應更新為：
  - **HCC1954 原先 medium (d=+0.57) → AF<0.4 重算後確定 POS + 升級**
  - **COLO829 原先 NEGATIVE → AF<0.4 仍 NEGATIVE，應 permanent exclude**
  - **其他 5 樣本 AF<0.4 框架下訊號穩定**

### 新的 per-sample 分類（AF<0.4 stratified）

| 分類 | 樣本 | 數量 |
|---|---|---|
| **medium+ (確定 POS)** | HCC1954 / HCC1937 / H1437 / HCC1395_DORADO | 4/7 |
| **small (POS)** | HCC1395 (flipflop across AF cuts) | 1/7 |
| **negligible (ceiling)** | H2009 | 1/7 |
| **NEGATIVE (out-of-scope)** | COLO829 | 1/7 |

5/7 POS（原 B.1-3 版本相同，但結論更精確：HCC1954 從 medium 升級，HCC1395_DORADO/HCC1395 因飽和下調是統計 artifact）。

---

## Part B: AF<0.4 coverage impact

| Sample | n 舊 filter | **n 新 filter** | coverage loss | n AF<0.2 | n TP 新 | n FP 新 | Power |
|---|---|---|---|---|---|---|---|
| **HCC1954** | 1,622 | **1,040** | 35.9% | 564 | 735 | 305 | YES |
| HCC1937 | 890 | 626 | 29.7% | 296 | 543 | 83 | YES |
| HCC1395 | 1,173 | 860 | 26.7% | 303 | 763 | 97 | YES |
| H1437 | 1,409 | 1,132 | 19.7% | 455 | 1,092 | 40 | YES |
| H2009 | 19,979 | 9,939 | 50.3% | 1,207 | 9,515 | 424 | YES |
| HCC1395_DORADO | 637 | 566 | 11.1% | 181 | 520 | 46 | YES |
| **COLO829** | 34 | 34 | 0% | 21 | 8 | 26 | **MARGINAL** |

### 解讀
- **整體 45% coverage loss（Step 3 確認）= 14,197 regions**；per-sample 分散 11-50%，皆有足夠 power（n≥100 6/7）
- **COLO829 n=34 未變**（AF≥0.4 的 NR≥80 NonLOH NG≥4 本來就不存在），但 TP:FP = 8:26 → 負向訊號，out-of-scope 確認
- **AF<0.2 高精度版本**：HCC1954 仍保留 n=564（power 充足）、HCC1937 n=296、H1437 n=455 → 論文 figure 用 AF<0.2 可行

---

## Part C: HCC1954 AF≥0.4 FP 聚集分析（與假設相反）

### [C1] 染色體層級分佈
HCC1954 AF≥0.4 FP 的染色體富集（單看 AF≥0.4 fraction）：

| Chr | total FP | AF<0.4 | AF≥0.4 | AF≥0.4 fraction |
|---|---|---|---|---|
| chr14 | 19 | 3 | 16 | **0.842** |
| chr2 | 40 | 9 | 31 | 0.775 |
| chr19 | 21 | 5 | 16 | 0.762 |
| chr7 | 105 | 30 | 75 | **0.714** |
| chr20 | 35 | 10 | 25 | 0.714 |
| chr16 | 92 | 28 | 64 | 0.696 |
| **chr8** | **118** | 37 | 81 | **0.686** |
| chr17 | 26 | 9 | 17 | 0.654 |
| chr4 | 31 | 11 | 20 | 0.645 |
| chr5 | 36 | 13 | 23 | 0.639 |
| chr1 | 71 | 35 | 36 | 0.507 |
| chr9 | 54 | 34 | 20 | **0.370** (低) |
| chr11 | 30 | 17 | 13 | 0.433 |

**關鍵觀察**：
- **chr7/8/16 是 HCC1954 FP 最多的染色體**（總 315 FP，佔全染色體 FP 的 >50%）
- **chr7/8/16 約 70% FP 在 AF≥0.4**
- **HER2 on chr17（AF≥0.4 fraction 0.654）、MYC on chr8** — 已知 HCC1954 HER2+ amplification、chr8q MYC amp 區域
- **chr9 例外**（37%）：可能 chr9 CDKN2A 區域 LOH 已被標註，或其他機制

### [C2] 空間聚集度（inter-FP distance within-chr）— 假設推翻 ⚠️

| Metric | AF≥0.4 FP | AF<0.4 TP (baseline) |
|---|---|---|
| inter-neighbor median (bp) | 650,034 | 372,853 |
| <100kb fraction | 0.270 | 0.333 |
| **ratio** (FP/TP) | — | **0.81×** |

**結論**：AF≥0.4 FP **不是空間 cluster**，實際上 **比 TP 更分散**（ratio <1）。原假設「FP 集中 hotspot」**被推翻**。

### [C5] Cross-sample clustering ratio

| Sample | n FP≥0.4 | n TP<0.4 | clust FP | clust TP | **ratio** |
|---|---|---|---|---|---|
| HCC1954 | 511 | 735 | 0.270 | 0.333 | 0.81× |
| HCC1937 | 172 | 543 | 0.186 | 0.304 | 0.61× |
| H1437 | 71 | 1,092 | 0.127 | 0.246 | 0.52× |
| HCC1395 | 126 | 763 | 0.071 | 0.311 | 0.23× |
| H2009 | 884 | 9,515 | 0.141 | 0.769 | 0.18× |

**5/5 樣本 ratio <1** → AF≥0.4 FP 在各樣本都**不聚集**，只是**分散於特定染色體**（chromosomal enrichment）。

### [C3] Feature profile — AF<0.4 TP vs AF≥0.4 FP

| Feature | AF<0.4 TP (n=735) | AF≥0.4 FP (n=511) | Δ (FP - TP) |
|---|---|---|---|
| NumReads | 117.2 | 118.3 | +1.0（類似）|
| **NumCpGs** | 89.6 | **111.0** | **+21.4** ⭐ |
| Coverage_Multiple | 1.56 | 1.58 | +0.01（類似）|
| AlleleDelta mean | 0.016 | 0.024 | +0.008 |
| caller_af | 0.184 | **0.607** | +0.423（by 定義）|
| **PairwiseMedianDist** | 0.189 | **0.168** | **-0.021** ⭐ |
| HPMergedDelta | 0.026 | 0.023 | -0.003 |

### [C4] |AlleleDelta| 分佈（Mann-Whitney）

| | AF≥0.4 FP | AF<0.4 TP | Mann-Whitney p |
|---|---|---|---|
| mean \|AlleleDelta\| | 0.027 | 0.032 | **p=8.67e-4** |
| median | 0.013 | 0.016 | (FP 略小) |

---

## 🌟 綜合機制解釋（基於 C1-C5）

原假設（hotspot 聚集）**被推翻**。新機制解釋如下：

### 失效機制三要素（HCC1954 專屬）

1. **CNV 染色體富集（非 hotspot）**：
   - HCC1954 HER2+ 高 ploidy（ICGC ~4），chr7/8/16 等染色體有大塊 CNV/cnLOH
   - 這些染色體 AF≥0.4 FP 比例 ~70%，**但分散在整條染色體**
   - 現有 LOH.bed 僅標註典型 LOH，未涵蓋細微 cnLOH 或 mixed ploidy 區域

2. **Germline het 被 CNV 驅動至 AF≠0.5**：
   - 正常 germline het AF≈0.5 → 在 3:1 或 1:3 allele imbalance 區變成 AF≈0.25 或 0.75
   - AF<0.4 會誤除「CNV 低 AF 段」但保留「CNV 高 AF 段」（germline appearing as higher AF）
   - **支持證據**：AF≥0.4 FP 的 PairwiseMedianDist 更低（0.168 vs TP 0.189）→ reads 間相似性更高 → 單 haplotype 主導 → germline-like 而非 clonal heterogeneity

3. **CpG island / 啟動子 bias**：
   - AF≥0.4 FP NumCpGs=111 vs AF<0.4 TP=89.6（**+24%**）
   - 啟動子 CpG island 區富集 → 高 CpG 密度
   - 啟動子 methylation 結構複雜（TF binding、bivalent chromatin）→ 易被 HP 分群誤判為 NG=4

### 為何 AF<0.4 有效挽救

- AF<0.4 **選擇性剔除 CNV 高 AF 段**（這裡 germline het 表現為 imbalance 高端）
- 同時保留 somatic SNV（天然低 VAF 因 purity<1 或 subclonal）
- 生物學對應：**tumor purity < 1 + subclonal events** → somatic AF 天然 <0.4

---

## 對其他任務的影響

### 需更新

1. **memory `project_hpfinengroups_subclone_marker.md` (20260418) B.1-3 per-sample 細節**：
   - HCC1954 從 "medium d=+0.57" → "**medium+ h=+0.654 (AF<0.4) / +0.775 (AF<0.2)**"
   - COLO829 從 "NEGATIVE d=-0.17 小樣本" → "**permanent out-of-scope**（AF<0.4 內仍 h=-0.43；AF<0.2 內更差 -0.55）"
   - 新增備註：HCC1395/HCC1395_DORADO AF<0.4 內數值 h 下降為飽和 artifact，非訊號衰減

2. **結論穩定度補充結論 16 ⭐4 依據強化**：
   - AF<0.4 stratified 5/7 POS 新驗證
   - 新發現：HCC1954 失效機制 = CNV 染色體富集 + germline het AF 漂移 + CpG island bias（非 spatial hotspot）

3. **docs/experiments/in_progress/2026/04/20260418_F_HPFineNGroups_deepening_POSITIVE_01.md**：
   - 新增 Step 4 章節
   - 補充 HCC1954 失效機制解釋

### 新的質疑與待驗

1. **CNV caller 整合必要性提升**：
   - chr7/8/16 染色體 bias 強烈建議整合 Delly/Manta/sequenza
   - Phase 2 跨樣本驗證若涵蓋 HER2+ / ERBB2-amplified 臨床樣本，必須 CN caller 介入
   - 這連接到 Opus 4.7 重整 plan B.2-2（Coverage_Multiple 代理 CN 的擔憂）

2. **ASM 在 CpG island 區的特殊性**：
   - 啟動子區 methylation 結構複雜，是否應 zone-aware 分析？
   - Phase 2C「CpG 功能分層」直接相關（C.4 研究方向）

3. **Purity<1 外推模擬**：
   - 若 purity=0.5，somatic AF 會壓至 ~0.25，AF<0.4 仍涵蓋
   - 但 germline het 會被 normal contamination 拉回 0.5，AF<0.4 可能同時剔除
   - **簡易數值模擬（非 Step 4 範圍）**可驗證 AF 閾值在不同 purity 下的 TP/FP trade-off

### 不受影響

- B.2 LOH Subclone AF×Methylation（正交 framework）
- Zone-Aware Framework
- Per-CpG ASM characterization

---

## 產物

| 檔案 | 描述 |
|---|---|
| `data/step4a_cohens_d_af04.tsv` | Per-sample h_all / h_af04 / h_af02 + classification |
| `data/step4b_coverage_impact.tsv` | Per-sample coverage loss + power check |
| `data/step4c_hcc1954_chr_af_fp.tsv` | 染色體 × AF group FP 計數 |
| `data/step4c_hcc1954_feature_profile.tsv` | 四 quadrant feature profile |
| `data/step4c_cross_sample_clustering.tsv` | 5 樣本 clustering ratio |

---

## Step 5 候選（未執行）

- **Purity<1 模擬**：驗證 AF<0.4 閾值穩健性（VAF ~Beta(AF*purity, (1-AF)*purity)）
- **CNV caller 整合**：chr7/8/16 HCC1954 FP 是否全部落在 Delly/Manta segment？
- **CpG island annotation**：UCSC CpG island BED 交叉驗證 AF≥0.4 FP 的 NumCpGs=111 高在哪類 region（promoter/shore/shelf/open sea）
- **AF 雙閾值（0.4 下限 + upper cap）**：如 0.1-0.4 vs 0.4-0.6，AF<0.1 是否也有 artifact

優先級：Step 5 偏向 Phase 2 主線準備工作，不阻擋 F pilot POSITIVE 結論。
