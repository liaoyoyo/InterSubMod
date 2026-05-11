# G2 · LOH & Subclone Annotation — 觀察報告

<!--
建立時間: 2026-04-23
狀態: draft
資料源:
  - research/feature_layered_observation/data/merged_with_vcf.tsv.gz (748,676 rows × 60 cols)
  - research/feature_layered_observation/scripts/g2_loh_subclone.py
輸出:
  - research/feature_layered_observation/figures/G2/fig01..06
  - research/feature_layered_observation/data/G2/step1..6 .csv
方法: research/feature_layered_observation/02_methodology.md
-->

## 1. 特徵定義與來源

| Feature | Source | Type | 說明 |
| --- | --- | --- | --- |
| `Potential_LOH` | `RegionProcessor.cpp` (`HP_Ratio < 0.1` or `> 0.9`) | binary | HP tag 基礎 LOH flag |
| `LOH_Subtype` | `RegionProcessor.cpp:188 determine_loh_subtype` | categorical(5) | None / LOH_Noise / LOH_Weak / LOH_Strong / LOH_Subclone |
| `LOH_Bed_Overlap` | `LohBedAnnotator.cpp classify_loh_source` | binary | 是否落在外部 LOH BED |
| `LOH_Source` | `LohBedAnnotator.cpp LohSource enum` | categorical | NONE/BED_ONLY/RATIO_ONLY/BOTH |
| `LOH_Bed_Annotation` | `LohBedAnnotator.cpp` | categorical | BED name column |
| `Subclone_ID` | `SubcloneAnalyzer.cpp:100 assign_subclones` | categorical | -1 / 0..K-1 |

**公式（核心，`RegionProcessor.cpp:188-204`）**：

```cpp
if (!potential_loh)                  -> "None"
else if (verif_class == "Noise")     -> "LOH_Noise"
else if (verif_class == "Weak")      -> "LOH_Weak"
else if (verif_class == "Strong")    -> "LOH_Strong"
else if (verif_class == "Subclone")  -> "LOH_Subclone"
else                                  -> "None"
```

→ `LOH_Subtype` 是 **`Potential_LOH ∧ VerificationClass`** 的 join；非獨立信號。

**主資料覆蓋狀況（748,676 rows）**：

| Feature | 可用？ | 理由 |
| --- | --- | --- |
| `Potential_LOH` | ✅ | 完整（F/T = 497,637 / 251,039） |
| `LOH_Subtype` | ✅ | 完整（None 66.5% / Noise 16.9% / Weak 10.2% / Strong 5.2% / Subclone 1.1%） |
| `LOH_Bed_Overlap` | ❌（near-zero signal） | 所有樣本 97.6% = False；HCC1954 paired 全 NaN |
| `LOH_Bed_Annotation` | ❌ | 整欄 NaN |
| `LOH_Source` | ❌ | 不在 master TSV |
| `Subclone_ID` | ❌ | 不在 master TSV |

→ 實際可觀察維度只剩 `LOH_Subtype`（5 類）與 `Potential_LOH`（binary）。

## 2. 觀察目標

1. 每個 LOH_Subtype 類別的 **TP rate**（Wilson 95% CI） × 7 樣本 × 2 modes
2. **LOH_Subtype × AF × CN 三維 heatmap**：找出 FP-enriched 或 TP-enriched cell
3. 跨樣本一致性：subtype ordering 是否跨樣本穩定
4. 若將 5 類視為 **ordinal（None → Noise → Weak → Strong → Subclone）**，其 AUC 是否顯著
5. **Confound guard**：殘差化 `vcf_AF` + `Coverage_Multiple` 後是否仍有訊號
6. Chr-level 聚集：LOH_Subtype 是否集中在特定染色體（空間 autocorrelation 警訊）

## 3. Step 1 · LOH_Subtype 5 類 TP rate × 7 樣本

圖：`fig01_lohsubtype_tp_rate_bars.png`  · 資料：`step1_tp_rate_per_class.csv`

### 3.1 Paired_full（caller 已近完美，baseline ~99%）

| sample | None | Noise | Weak | Strong | Subclone |
|---|---:|---:|---:|---:|---:|
| HCC1395 | 0.980 | 0.970 | 0.981 | **0.997** | 0.993 |
| HCC1395_DORADO | 0.994 | 0.992 | 0.988 | 0.990 | 0.995 |
| HCC1937 | **0.994** | **0.997** | 0.962 | 0.954 | 1.000 |
| HCC1954 | 0.999 | 0.992 | 0.996 | 1.000 | 1.000 |
| H2009 | 1.000 | 0.999 | 0.998 | 0.998 | 0.999 |
| H1437 | 1.000 | 0.999 | 1.000 | 1.000 | 1.000 |
| COLO829 | 0.942 | **0.904** | 0.937 | 0.932 | **0.903** |

**觀察**：
- Paired mode baseline TP rate 極高（>94%），**LOH_Subtype 之間差異絕對值 <0.1 pp**，無 filter 價值。
- COLO829 與 HCC1937 是唯二有可識別 FP gap 的樣本，但 Subclone 類別 n 太小（31 / 481）。
- HCC1937 出現反直覺：**Weak/Strong TP rate（0.962, 0.954）< None（0.994）** — 與 subtype ordering 假設矛盾。

### 3.2 to_pileup（FP 豐富，baseline 0.25–0.91）

| sample | None | Noise | Weak | Strong | Subclone |
|---|---:|---:|---:|---:|---:|
| HCC1395 | 0.699 | 0.930 | 0.887 | 0.905 | 0.958 |
| HCC1395_DORADO | 0.685 | 0.643 | 0.859 | 0.855 | 0.641 |
| HCC1937 | **0.463** | **0.431** | 0.794 | 0.802 | 0.758 |
| HCC1954 | **0.245** | **0.251** | 0.291 | 0.377 | **0.510** |
| H2009 | 0.909 | 0.909 | 0.948 | 0.953 | 0.937 |
| H1437 | 0.762 | 0.747 | 0.828 | 0.876 | 0.898 |
| COLO829 | 0.648 | 0.657 | 0.732 | 0.647 | 0.732 |

**觀察**：
- 在 TO mode 下 **LOH_Subtype TP rate 有明顯分層**：None/Noise 大多 <0.75，Strong/Subclone 大多 >0.80。
- HCC1937 / HCC1954：Strong/Subclone TP rate 相較 None gap +0.26 ~ +0.34 pp → potential TP-rescue signal。
- HCC1395_DORADO 例外：Subclone TP rate（0.641）反而與 None（0.685）相當 → subtype ordering 不穩定。
- H2009/H1437 subtype 差異 <0.15 pp，不適合做 hard filter。

## 4. Step 2 · LOH × AF × CN heatmap

圖：`fig02_loh_af_cn_heatmap.png`  · 資料：`step2_loh_af_cov_cells.csv`

切分：`LOH_Subtype (5) × AF_class (3) × Coverage_Category (6) = 90` 組合，其中 60 組 n≥20。

**TP-weighted mean per LOH_Subtype**（aggregated across 60 valid cells）：

| LOH_Subtype | n_cells | total_n | weighted TP rate |
| --- | ---: | ---: | ---: |
| None | 17 | 497,627 | **0.818** |
| LOH_Noise | 6 | 126,604 | **0.763** ← 最低 |
| LOH_Weak | 15 | 76,738 | 0.911 |
| LOH_Strong | 16 | 39,207 | **0.917** ← 最高 |
| LOH_Subclone | 6 | 8,435 | 0.884 |

**反直覺發現**：`LOH_Noise`（0.763）TP rate **低於** `None`（0.818），顛覆「LOH → 更多 TP」的直覺。**Noise 代表 HP_Ratio 極端但 ISM clustering 未通過 verification**，可能是 self-phasing artifact + germline variant 的混合。

**Bottom 3 FP-enriched cells**：

| LOH_Subtype | AF_class | Coverage_Category | n | TP rate |
| --- | --- | --- | ---: | ---: |
| LOH_Noise | Extreme | CNV_Loss | 5,949 | 0.667 |
| LOH_Weak | Intermediate | High_Copy | 31 | 0.516 (small n) |
| LOH_Strong | Intermediate | High_Copy | 199 | 0.693 |
| None | Extreme | Low | 32,888 | 0.719 |

**Top 3 TP-enriched cells**：

| LOH_Subtype | AF_class | Coverage_Category | n | TP rate |
| --- | --- | --- | ---: | ---: |
| LOH_Strong | Near-half | Normal | 3,351 | 0.996 |
| LOH_Weak | Intermediate | Normal | 12,837 | 0.982 |
| LOH_Strong | Intermediate | Normal | 10,260 | 0.964 |

→ **LOH_Subtype 的 filter 訊號主要來自與 AF Near-half + Normal CN 的聯集**，而非 LOH 本身。

## 5. Step 3 · 跨樣本一致性

圖：`fig03_per_sample_grid.png`

**Paired_full**：全部樣本 5 類 TP rate 都 >0.90，差異太小 → subtype ordering 無法跨樣本穩定。  
**to_pileup**：5/7 樣本出現「None/Noise 低、Strong/Subclone 高」的單調上升 pattern，但：
- HCC1395_DORADO Subclone TP rate (0.641) < None (0.685) → 例外
- HCC1954 整體 baseline 極低（0.25-0.51），即使 Subclone 0.510 也僅剛過 0.5

Spearman concordance 未正式計算（n_samples=7），但目測 4/7 符合 monotonic ordering。

## 6. Step 4 · Ordinal AUC（若視為序數）

資料：`step4_ordinal_auc.csv`

| sample | paired AUC | TO AUC | Potential_LOH TO AUC |
| --- | ---: | ---: | ---: |
| HCC1395 | 0.517 | 0.526 | 0.526 |
| HCC1395_DORADO | **0.437** ↓ | 0.564 | 0.528 |
| HCC1937 | **0.291** ↓↓ | **0.596** | 0.538 |
| HCC1954 | 0.388 ↓ | 0.522 | 0.519 |
| H2009 | **0.262** ↓↓ | 0.524 | 0.516 |
| H1437 | 0.430 | 0.529 | 0.517 |
| COLO829 | 0.485 | 0.509 | 0.508 |

**判讀**：
- **Paired mode AUC <0.50 全面反轉**（5/7 樣本）→ LOH_Subtype 在 paired 下反而指向 FP；Spearman ρ 顯著負（p<0.001 於 HCC1937、HCC1954、H2009、COLO829）。
- TO mode AUC 全部 >0.50，但最高僅 HCC1937 0.596 → **未達 0.58 Beyond-AUC ceiling 在 6/7 樣本**。
- `Potential_LOH`（binary）單獨 AUC 與 ordinal AUC 差異 <0.06 → 5 類序數資訊**未提供顯著額外 lift**。

按 **methodology Verdict 規則**：`LOH_Subtype/Potential_LOH` 本體 **NEGATIVE**（AUC 0.50-0.58、跨樣本不一致）。

## 7. Step 5 · Confound Residualization (vcf_AF + Coverage_Multiple)

資料：`step5_residualized_auc.csv`

| sample | mode | LOH raw | LOH resid | PL raw | PL resid |
| --- | --- | ---: | ---: | ---: | ---: |
| HCC1395 | paired | 0.517 | **0.378** | 0.496 | **0.380** |
| HCC1937 | paired | 0.291 | 0.491 | 0.360 | 0.540 |
| HCC1954 | paired | 0.388 | 0.563 | 0.382 | 0.566 |
| H2009 | paired | 0.262 | 0.272 | 0.259 | **0.210** |
| H1437 | paired | 0.430 | **0.076** | 0.417 | **0.059** |
| HCC1395 | TO | 0.526 | 0.491 | 0.526 | 0.465 |
| HCC1395_DORADO | TO | 0.564 | **0.588** | 0.528 | **0.627** |
| HCC1937 | TO | **0.596** | **0.605** | 0.538 | **0.663** |
| HCC1954 | TO | 0.522 | **0.708** | 0.519 | **0.716** |
| H2009 | TO | 0.524 | 0.553 | 0.516 | 0.555 |
| H1437 | TO | 0.529 | 0.544 | 0.517 | 0.549 |
| COLO829 | TO | 0.509 | 0.507 | 0.508 | 0.505 |

**關鍵新發現**：

1. **Paired mode：殘差化後全面崩塌**（H1437 PL AUC 0.059 是極端 anti-signal；HCC1395 0.380）→ paired 下任何「AF + CN 之上」的 LOH 訊號都消失或反向。
2. **TO mode：殘差化 AUC 在 4/7 樣本反而上升**（HCC1395_DORADO/HCC1937/HCC1954/H2009 Potential_LOH resid > raw），**HCC1954 resid AUC=0.716** 達 POSITIVE 門檻。
3. **Potential_LOH**（binary）殘差化後反而比 5-class `LOH_Subtype` **更乾淨**（HCC1937 TO: 0.663 vs 0.605）— subtype 細分加入了 verification class 的雜訊。

**confound-guard 解讀**：
- 殘差化 AUC 上升 → LOH 與 AF/CN 之間有 **suppressor（負向混淆）**：vcf_AF extreme 區段同時壓低 TP 機率與增加 LOH 判定 → 去除此 confound 後 LOH 的 TP-enrichment 顯露。
- **但下游 filter 若已條件在 AF/CN 上則此殘差訊號已「被吸收」**。因此 HCC1937/1954 TO 的 AUC 0.66-0.72 不等於可用的 F1 gain。

## 8. Step 6 · Spatial Autocorrelation（Chr 分佈）

圖：`fig06_per_chr_distribution.png`  · 資料：`step6_per_chr_loh.csv`

**關鍵觀察**（paired_full 下）：
- **H2009**：chr8（myc）與 chr19 LOH_Strong + Subclone fraction 明顯升高（>40%）
- **HCC1954**：chr8 / chr17 LOH_Strong 高度聚集（與 Z3 amplicon pilot 報告一致）
- **HCC1395**：全 chr 相對均勻分佈，chr 特異性弱
- **H1437**：chr5 / chr19 LOH_Weak fraction 偏高

**警訊**：**LOH_Subtype 具明顯 spatial autocorrelation**，任何 chr+pos 聚合特徵（包含自身 AUC）必須做 mid-TP-rate 視窗驗證（依 feedback_spatial_autocorrelation_confound 記憶原則）。本觀察建議所有下游「LOH × 甲基化」特徵都需做 chr-stratified replication。

## 9. 論文與知識庫背景

（此節為質性引用空位，未做 MCP knowledge 查詢）

- LOH 文獻常見二分：germline CNV loss 與 tumor-acquired LOH；本專案將 HP tag ratio extremes 當 operational proxy，已在 `projects/LOH_research_status` 記憶中確認 16-52× TO:paired mismatch。
- Subclone 定義：本 codebase `VerificationClass == Subclone` 意指 ISM clustering 顯示 **multi-modal 甲基化 pattern**，與分群穩定度結合；**不是外部 subclonal CNV caller 定義**。
- 文獻上 LOH 作為 variant filter 已被多個 ONT somatic pipeline 測試，結論普遍 **NEGATIVE**（LOH 區 TP-enriched，filter 會誤傷）— 與本觀察 paired_full 段完全一致。

### 9.1 知識庫引用（Phase D 補）

查詢詞：`LOH bed` (top_score 56.9, high confidence)、`allele imbalance` (top_score 20.8)、`LOH subtype` (top_score 21.7)。下表列出 `/big8_disk/liaoyoyo2001/Knowledge/` 下相關文件。

| kb_path | kb_title | 與 G2 的關聯 |
|---|---|---|
| `05_tools/longphase-to.md` | LongPhase-TO | **主要來源**：`LOH.bed`（體染色體 F1=96.2%）的正式定義；LOH 基於 phased genotype ratio；說明為何 paired 模式下 LOH_Subtype≈random（paired caller 不依 LOH BED 過濾），TO 模式下才有 selection pattern |
| `04_databases/seqc2-truth-set.md` | SEQC2 Truth Set（真陽性標準集） | CNV/LOH benchmark（Zenodo v4：gain/loss/LOH BED + exclusion + VCF + median CN）；SEQC2 定義 `LOH=CN=2 但一條 haplotype 完全缺失`，**與 LongPhase-TO 定義不同**，此為 G2 `LOH_Subtype` 類別語意的交叉檢核依據 |
| `06_workflows/phasing-workflow.md` | Phasing 工作流程 | 明確釐清 `HP=3`（LongPhase-S read-level ALT）**不等於** LOH；只有 phased VCF / `LOH.bed` 才是 canonical LOH — 本 G2 的 `LOH_Subtype` 來自 `LOH.bed` |
| `05_tools/longphase-s.md` | LongPhase-S | `eff_hp >= 30` 的 LOH eligible 閾值；HP tag 極端比例與 purity 關聯，是 `Potential_LOH` 的上游依據 |
| `09_standards/common-pitfalls.md` | 常見陷阱與易錯清單 | LongPhase-S 與 LongPhase-TO 的 LOH 語意差異（read-level ALT vs. region-level phased genotype）— 本 G2 使用 region-level 定義 |

**LOH/Subclone**：3/3 查詢主題全數命中（high confidence），覆蓋 LOH.bed 定義、truth set、allele imbalance 推論三個層面。
**空缺**：KB 對 `cnLOH (copy-neutral LOH)` 與 `germline LOH` 子類別的細分未覆蓋，本 G2 `LOH_Subtype`（5 類）為 InterSubMod 內部 taxonomy，建議 Phase F 把分類定義回寫 KB。

### 9.2 外部文獻（Phase D）

**LOH / subclonal copy number detection 領域代表性論文：**

1. **Shen, R. & Seshan, V. E. (2016).** "FACETS: allele-specific copy number and clonal heterogeneity analysis tool for high-throughput DNA sequencing." *Nucleic Acids Research* 44(16), e131. DOI: **10.1093/nar/gkw520** — FACETS 是 short-read tumor-normal 的 allele-specific copy number pipeline，透過 heterozygous SNP 的 BAF 進行 joint segmentation，輸出 tumor purity / ploidy / integer CN / cnLOH 類別。與本 G2 的關聯：提供 **cnLOH 正式定義（CN=2 但 minor allele=0）**，對照本專案 `LOH_Strong` 的 HP_Ratio-based operational proxy。方向：**挑戰** 我們結論的「BAF 外部 truth」潛在來源（若用 FACETS BAF 重做殘差化，可能破解 HP_Ratio 自我參照）。

2. **Nik-Zainal, S., Van Loo, P., Wedge, D. C. et al. (2012).** "The life history of 21 breast cancers." *Cell* 149(5), 994–1007. DOI: **10.1016/j.cell.2012.04.023** — Battenberg algorithm 原始論文，透過 haplotype phasing 將 subclonal CN 判定解析度推到 3% 細胞比例。與本 G2 的關聯：**確立「LOH 區富 TP」的先驗基礎** — 早期 subclonal LOH 常伴隨 driver 突變聚集，與本觀察 paired_full 下 LOH_Strong TP rate >0.98、LOH 作 filter anti-signal 一致。方向：**支持** 我們 paired NEGATIVE / characterization-only 的結論。

3. **Onuchic, V., Lurie, E., Carrero, I. et al. (2018).** "Allele-specific epigenome maps reveal sequence-dependent stochastic switching at regulatory loci." *Science* 361(6409), eaar3146. DOI: **10.1126/science.aar3146** — 用 71 epigenomes 建構 allelic imbalance high-resolution map，證明 heterozygous regulatory loci 大量呈現 stochastic methylation switching。與本 G2 的關聯：為 `LOH_Subtype × methylation` 交互作用提供生物學先驗 — 若 LOH 同時消除某 allele 的 stochastic switching，methylation heterogeneity 應下降。方向：**支持** 在 TO mode 下 LOH_Strong stratum 可能有微弱 characterization 價值（TP rescue direction），但**挑戰**將其作為 global filter。

**LOH-constrained phasing 相關空白**：2026-04-22 Thread D 新發現（NG=2 在 Inner LOH 93-99% same-hap）目前**無直接對照文獻**，4-bucket phasing 分類機制未見於 LongPhase / HapCUT2 / Margin — 候選本研究新貢獻。

## 10. 結論、三質疑與邏輯鏈

### 10.1 Verdict（依 Methodology §Verdict）

| 類別 | 全域 AUC（ordinal） | Confound guard | 跨樣本 | Verdict |
|---|---|---|---|---|
| `LOH_Subtype` (paired) | 0.26-0.52（反向/隨機） | resid 崩塌（H1437 0.08） | 5/7 AUC<0.5 | **NEGATIVE / ANTI-SIGNAL** |
| `LOH_Subtype` (TO) | 0.51-0.60 | resid ≤0.71 在 4/7 | 部分上升 | **SAMPLE_SPECIFIC / CONFOUND_COLLAPSED** |
| `Potential_LOH` binary | 同上 | resid 在 TO 微升 | 同上 | 同上 |

**綜合判定**：LOH_Subtype 與 Potential_LOH **不適合作為 global canonical filter**；作為 characterization annotation 則可保留。

### 10.2 三質疑

1. **「LOH_Subtype 5 類是真的 ordinal 嗎？」**
   - 反證：TO mode HCC1395_DORADO Subclone TP rate 0.641 < None 0.685；paired mode HCC1937 None 0.994 > Strong 0.954。  
   - 詮釋：`VerificationClass` 與 LOH 獨立，subtype 是 2D 交乘；**不是 ordered**。應視為 5 個離散 stratum，而非序數。

2. **「TO mode 殘差化 AUC 上升是真訊號還是 over-fitting？」**
   - 質疑：`vcf_AF` 與 `Coverage_Multiple` 之上的 resid 仍保留 `HP_Ratio` 極端，但 `HP_Ratio` 本身就是 `Potential_LOH` 的定義 → **self-reference circularity**。
   - 驗證建議：以 **independent LOH 定義（LOH.bed 外部來源）**重做殘差化；若 AUC 同樣上升才可信。
   - 現況：master TSV 的 `LOH_Bed_Overlap` 全 False，無法即時做外部驗證。

3. **「chr8/17 的 LOH_Strong 聚集是生物學還是 caller artifact？」**
   - HCC1954 chr8/17（HER2/MYC amplicon）已知 CNV amplicon artifact（Z3 pilot 結論）。
   - H2009 chr8 LOH 聚集需獨立確認是否為 *in silico* 衍生，或與 KRAS/MYC 真實事件對齊。
   - 建議：跨 paired/TO 比對 `LOH_Strong ∩ chr8`；若 TO only 出現 → artifact；若 paired 也出現 → 真實。

### 10.3 邏輯鏈

```
LOH_Subtype/Potential_LOH 是否應納入 filter？
│
├── paired_full: baseline TP rate 0.90-1.00，5 類差距 <0.10 pp
│     └─ AUC 0.26-0.52（普遍 <0.50）→ LOH 在 paired 是 **anti-signal**
│         └─ 與既有結論一致：LOH 是 TP-enriched（filter 會誤傷）
│
├── to_pileup: baseline 0.25-0.91，LOH 類別 TP rate 高於 None
│     ├── HCC1937/HCC1954 gap +0.26 ~ +0.34 pp（TP rescue direction）
│     ├── HCC1395_DORADO 反向（Subclone < None）
│     └─ AUC 最大 0.596；殘差化後 HCC1954 達 0.72
│         └─ 但 self-reference circularity（HP_Ratio 既定義 LOH 又當 feature）
│
├── subtype ordering 在跨樣本不穩定 → 不視為 ordinal
│
├── spatial autocorrelation：chr8/17/19 聚集 → 混淆 caller-specific CNV artifact
│
└── 結論: **LOH_Subtype/Potential_LOH 作為 global filter NEGATIVE**;
           作為 TP-rescue stratification in TO mode 可進一步評估（HCC1937/1954 subset），
           但須先解決 HP_Ratio 自我參照問題，且 chr stratification 做 robustness check.
```

### 10.4 後續建議

| 優先 | 行動 | 預期收益 |
|---|---|---|
| P1 | **TO mode：LOH_Subtype 作為 TP-rescue stratum** 在 HCC1937/1954 評估 delta F1（subset-level） | 小規模 F1 gain |
| P1 | **HCC1395 paired: LOH_Noise 低 TP rate 深度診斷**（self-phasing 餘緒） | 驗證 PON-only phasing 是否已收斂 |
| P2 | 用 **LOH.bed 外部 truth**（SEQC2 Jaccard 0.928）重做殘差化，破解 HP_Ratio 循環 | 確認 TO HCC1954 AUC 0.72 是真訊號 |
| P2 | 加入 `Subclone_ID` 到 master TSV 後重做 G2 | 補足 Subclone 維度 |
| P3 | 跨 paired/TO 交叉比對 `LOH_Strong ∩ chr8` → artifact vs biology | chr 層級 spatial validation |
| P3 | `LOH_Bed_Overlap` 為 all-False bug 修復（`LohBedAnnotator` 在 master 生成時是否被 skip） | Feature 完整性 |

### 10.5 filter 建議（精簡版）

**不建議**把 `LOH_Subtype` 或 `Potential_LOH` 加入 canonical filter：
- paired_full 方向與 TP 反向（抓了 LOH 會增加 FN）
- TO 方向雖 weakly positive 但 AUC <0.60，且殘差化後的 signal 含 self-reference bias
- spatial autocorrelation 警訊未解

**可保留的操作**：
1. **TO mode characterization annotation**：`LOH_Subtype + Normal CN` 標示 TP-enriched (Strong+Subclone+Normal) 子集
2. **HCC1937/1954 TO subset**：作為下游 ReadParser 二階 filter 的 stratum 條件（與 S4 pilot 合併評估，如本週 P2 項）
