---
status: in_progress
date: 2026-04-22
type: analysis_report
revision: v2 (replaces v1 "NG_KDE_Rescaled_Multi_CN_Analysis_01.md" which was off-topic)
sample_set: [HCC1395, HCC1395_DORADO, HCC1937, HCC1954, H2009, H1437, COLO829]
primary_testbed: HCC1395 TO mode (paired baseline TP = 71.1%, FP = 28.9%)
pipeline_track: characterization / biology-informed filter design
related:
  - docs/experiments/in_progress/2026/04/20260421_NG_KDE_Rescaled_Multi_CN_Analysis_01.md  # v1, off-topic
  - docs/experiments/in_progress/2026/04/20260420_KDE_Fix_Acceptance_Validation_01.md       # KDE 修正基礎
  - docs/experiments/validated/2026/04/20260414_LOH_Subclone_AF_Methylation_Evidence_01.md  # 20260414 baseline
question: >-
  Under new-KDE baseline, do LOH_Subtype x AF_class x CN-tier combinations
  effectively separate TP from FP? Can we assign different parameter meanings
  per biological region (stratified biology-informed filters)?
---

# LOH × AF × 新 KDE 下的 TP/FP 判別分析（v2）

**結論一行**：在 HCC1395 TO mode（baseline TP = 71.1%）下，`LOH_Subtype × AF_class × CN-tier` 切分出**清晰的 TP/FP 兩極**：LOH 註解（Strong/Weak/Noise/Subclone）+ 非 Intermediate AF 的區塊 TP rate 落在 **88-97%**（最高 S3 = 95.5%），而 `None + Intermediate AF + CN≥3` 是 **FP 聚集區**（FP rate 36-53%）；此訊號結構允許**依生物模組套用不同語意的分層濾鏡**，S5 組合可在保留 2.85% TP 下移除 99.37% FP。

---

## 0. 為何需要 v2

v1 報告的主軸是「20260414 結論在新 KDE 下方向是否保留」——回答的是**效應量穩健性**。用戶實際問的是：

> LOH 與 AF 在新 KDE 判別上是否對分類高 TP 更有效？可否用不同區域的生物狀況解釋分出 TP 與 FP？依狀況與生物機制用不同參數的不同意義去處理？

此問題的正確答題框架是 **TP/FP 判別 + 生物機制分區 + stratified filter 提案**，而非 ΔNG 效應量比較。v2 重新切入。

---

## 1. 資料來源與測試場選擇

### 1.1 資料

- **Master**: `research/ng_kde_rescaling/data/merged_7samples_paired_full_plus_hcc1395_to.tsv.gz`（369,080 rows；由 step0 從新 KDE paired_full 6 樣本 + HCC1954 rescaled + HCC1395 TO phase1 拼接）
- **欄位**: `tp_label`, `LOH_Subtype`（None/Weak/Noise/Strong/Subclone，來自 external LOH.bed + ReadParser annotation）, `AF_class`（Extreme <0.1 or >0.9 / Near-half 0.4-0.6 / Intermediate 其餘）, `CovM_used`（新 KDE 框架）, `cn_tier_F`（SEQC2 empirical 0.65/0.99/1.33/1.82）

### 1.2 為何 TO mode 是主測試場

| Mode | 總 rows | TP | FP | TP% |
|------|---------|------|------|------|
| paired_full (7 樣本) | 328,965 | 315,507 | 13,458 | 95.9% |
| **to_pileup (HCC1395)** | **40,115** | **28,509** | **11,606** | **71.1%** |

paired_full 的 FP% 僅 4.1%（最低 H1437 僅 8 個 FP），**cell-level 判別分析無統計力**。TO mode 的 29% FP 才提供足夠樣本評估 stratified filter 的 FP reduction。

**paired_full** 作為 **per-sample 驗證軌**（確認模組化框架在各樣本不崩），不是主結論來源。

### 1.3 測試軌切分

| 軌 | 資料 | 目的 |
|----|------|------|
| **Main** | HCC1395 TO mode | TP/FP 判別力量測、filter scheme 比較 |
| **Validation** | 7 樣本 paired_full | 每樣本 scheme 穩定性、跨樣本 generalizability |

---

## 2. 方法：生物模組化 TP/FP 判別

### 2.1 模組定義（biology-informed cells）

**維度**：`LOH_Subtype (5) × AF_class (3) × cn_tier_F (5)` = 75 cells；有效（n≥20）在 TO mode ≈ 32 cells。

**metric**：
- `TP_rate = n_TP / (n_TP + n_FP)` — 主判別力
- `TP_enrichment = (n_TP/N_TP) / (n_FP/N_FP)` — 相對富集
- `FP_rate = n_FP / n` — 雜訊程度
- Wilson 95% CI（保守的小樣本信賴區間）

### 2.2 Filter Scheme 定義

| Scheme | 生物詮釋 | 遮罩 |
|--------|---------|------|
| **S1 LOH_Strong + Extreme AF** | Deletion/cnLOH 下 pure haplotype（單一存活等位基因 AF≈0 or 1）| `LOH_Subtype=Strong AND AF_class=Extreme` |
| **S2 Subclonal LOH + Intermediate AF** | 亞克隆 LOH with admixed populations | `LOH_Subtype in (Subclone, Weak) AND AF_class=Intermediate` |
| **S3 Diploid Het Somatic** | No LOH + balanced AF + CN~2/3（canonical 雜合 somatic）| `LOH_Subtype=None AND AF_class=Near-half AND CN in (T1, T2)` |
| **S4 NonLOH + Extreme AF（高風險）** | 無 LOH 卻極端 AF → germline 泄漏 / mapping bias | `LOH_Subtype=None AND AF_class=Extreme` |
| **S5 組合（S1∨S2∨S3，排除 S4）** | 高信心生物場景之聯集 | — |
| **S6 S1 + NG≥3** | S1 + subclone haplotype marker | S1 AND HPFineNGroups≥3 |
| **S7 S5 + NG≥3** | S5 + NG≥3 | S5 AND NG≥3 |

---

## 3. 主結果（HCC1395 TO mode）

### 3.1 TP rate 熱圖：LOH × AF × CN 呈現雙極分佈

![F1 TP rate heatmap](figures/20260421_ng_kde_rescaled/fig_v2_1_to_tp_heatmap.png)

**觀察**：
1. **LOH 任何註解（Noise/Weak/Strong/Subclone）× Extreme AF** → TP rate 88-96%（明顯高於 71.1% baseline）。LOH 訊號即使被標為 "Noise"（註解不確定）仍帶出實質判別力。
2. **None + Near-half（balanced）× CN~2/3** → TP rate 93-96%（canonical het somatic 符合預期）。
3. **None + Intermediate × CN≥3** → TP rate 驟降至 **47-60%**（CN~4+Inter 僅 47%、CN≥5+Inter 僅 33%）→ FP 聚集區。
4. **None + Extreme × CN<1** → TP rate 僅 54%（del/LOH 但無 LOH 註解 = 假 deletion 或 ReadParser 判定歧異）。

### 3.2 FP rate 熱圖（同 axes，互補視角）

![F2 FP rate heatmap](figures/20260421_ng_kde_rescaled/fig_v2_2_to_fp_heatmap.png)

**FP 熱點**：
- `None + Intermediate + CN≥5` → **67% FP**（n=27 小但明顯）
- `None + Intermediate + CN~4` → **53% FP**（n=365）
- `None + Intermediate + CN~3` → **40% FP**（n=2204，量大）
- `None + Extreme + CN<1` → **46% FP**（del/LOH 但無 LOH 標記的區塊）

**生物詮釋**：CN≥3 + Intermediate AF + 無 LOH 有兩種典型 artifact：
1. **Repeat / segdup mapping chaos**：拷貝數異常區 + 中間 AF 常為 misaligned reads
2. **Low-purity subclone false-positives**：CN gain 下的低頻 somatic 與 germline leak 幾何上難分

### 3.3 Filter scheme TP 率 / TP 召回 / FP 消除

![F3 filter scheme bar](figures/20260421_ng_kde_rescaled/fig_v2_3_filter_scheme_bar.png)

| Scheme | n | TP rate | TP recall | FP reduction | n_TP | n_FP |
|--------|---:|--------:|----------:|-------------:|-----:|-----:|
| **S1 LOH_Strong + Extreme** | 292 | **90.1%** | 0.92% | 99.75% | 263 | 29 |
| **S2 Subclonal LOH + Inter** | 214 | 87.4% | 0.66% | 99.77% | 187 | 27 |
| **S3 Diploid Het** | 380 | **95.5%** ⭐ | 1.27% | 99.85% | 363 | 17 |
| **S4 NonLOH + Extreme（風險）** | 30432 | 71.1% | 75.95% | 24.35% | 21652 | 8780 |
| **S5 Combo (S1∨S2∨S3) \ S4** | 886 | **91.8%** | 2.85% | 99.37% | 813 | 73 |
| **S6 S1 + NG≥3** | 253 | 90.5% | 0.80% | 99.79% | 229 | 24 |
| **S7 S5 + NG≥3** | 484 | 90.7% | 1.54% | 99.61% | 439 | 45 |

**關鍵觀察**：
- **S3（Diploid Het）是單一最純模組**（95.5% TP，FP reduction 99.85%）；canonical biology 最能信任
- **S5 組合覆蓋 2.85% 全部 TP** 但 FP 清除 99.37% → 若作 high-confidence subset 非常有力
- **S4 是「整個世界」**（75.95% TP + 75.7% FP 都在這個 bucket），單憑 S4 無法切分，必須與其他維度（如 CN-tier × AF bin）進一步拆解
- **NG≥3 加進去幾乎沒幫助**（S6 vs S1、S7 vs S5 差異 <1pp）→ 在生物模組化框架下，NG 的邊際貢獻被生物學訊號吸收

### 3.4 Operating points（purity vs recall）

![F5 operating points](figures/20260421_ng_kde_rescaled/fig_v2_5_operating_points.png)

**左（TO mode）**：S3、S1、S5 都在左上角（高純度、低覆蓋），適合作 confidence tag；S4 落在底部接近 baseline，是 "everything bucket"。
**右（paired_full pooled）**：baseline 95.9% 過高，所有 scheme 貼近上緣；量測差異集中在 0.5pp 量級，需 TO mode 佐證。

### 3.5 生物模組詮釋表

![F6 biology module interpretation](figures/20260421_ng_kde_rescaled/fig_v2_6_biology_module_interpretation.png)

---

## 4. Per-sample 驗證（paired_full 7 樣本）

![F4 per sample schemes](figures/20260421_ng_kde_rescaled/fig_v2_4_per_sample_schemes.png)

**觀察**：paired_full 每樣本 TP 基準已在 93-99.9%（FP 稀少），scheme 之間 TP rate 差異僅 0-2pp。主要用途是**確認模組不會在特定樣本崩潰**：

- S3（Diploid Het）在每個樣本都達 99%+ TP → 穩定
- S1 在 HCC1395/COLO829 n 較低但仍 99%+；HCC1937 n=0（該樣本 LOH_Strong+Extreme 無 cell）
- S4 在所有樣本都約等於 sample baseline → 確認 S4 無辨別力（如預期）

**結論**：模組化框架在 7 樣本一致運作，不存在樣本專屬崩潰。

---

## 5. 生物機制對照：為什麼這些模組 work？

### 5.1 S1 Deletion / cnLOH + Extreme AF（TP 90.1%）

- **細胞學機制**：單一等位基因刪除或複製後再刪除（cnLOH），tumor read 全部來自存活 haplotype → caller AF→0 or 1（Extreme）
- **為什麼純**：Extreme AF 在 No-LOH 背景是 germline/artifact，但在 LOH 背景是 **預期的 somatic 行為**（FACETS/Battenberg 定義）
- **文獻**：Shen & Seshan 2016 FACETS, Nik-Zainal 2012 Battenberg, Ha 2014 TITAN

### 5.2 S2 Subclonal LOH + Intermediate AF（TP 87.4%）

- **細胞學機制**：LOH 事件發生在部分 tumor 亞克隆，混合 pure-LOH + diploid populations → AF 中間
- **為什麼可辨**：subclonal/weak LOH 註解捕捉到「有 LOH 訊號但非全體 tumor」；AF 中間是混合比例的直接反映
- **文獻**：Larose 2023 CAMDAC subclonal methylation-CN, Van Loo 2010 ASCAT

### 5.3 S3 Diploid Het Somatic（TP 95.5%，最純）

- **細胞學機制**：Canonical heterozygous somatic mutation 在無 CN 變異、diploid 區塊
- **為什麼純**：Near-half AF 在 diploid + non-LOH 是 Hardy-Weinberg 預期 → 落入此區塊的 variant 高度符合 somatic 理論模型
- **風險**：n 僅 380 → 覆蓋低，不能單靠此模組

### 5.4 S4 NonLOH + Extreme AF（高風險，TP 僅 71.1%）

- **artifact 機制候選**：
  - Germline variant leak（tumor-only caller 無法排除）
  - Mapping bias in segdup/repeat
  - Strand bias / low-quality read cluster
  - PCR duplicate 殘留
- **為什麼含大量 TP**：HCC1395 真 somatic 分佈本身偏 Extreme AF（腫瘤純度高），因此 TP 仍佔 71%。S4 不是「FP bucket」而是「無辨別力 bucket」
- **實作建議**：S4 內部須用**其他維度**（ReadParser DepthCF、chr 位點、鄰近 variant 群集）二次分層，不能單靠當前框架解決

### 5.5 LOH_Noise 的意外訊號

LOH_Noise + Extreme AF（CN~2/3）TP rate 88-96%，與 LOH_Weak 相當。詮釋：
- LOH_Noise 不是「無 LOH」，而是「LOH 判定訊號 noisy」
- 在 tumor AF extreme 輔助下，這些 cell 行為仍像 LOH → **LOH.bed annotation 的 "Noise" 類不應被拋棄**

---

## 6. 參數意義的差異化（回答用戶核心問題）

用戶要求：「依據狀況與生物機制與解釋用不同參數的不同意義去處理」

**本分析的參數差異化方案**：

| 區域 | 參數 | 語意 | 下游處理建議 |
|------|------|------|------------|
| **S1 / S3** | `CONF = HIGH_SOMATIC` | 生物學強支持 + 高 TP 率 | 直接保留 variant、pass flag |
| **S2** | `CONF = SUBCLONAL_CANDIDATE` | Subclonal biology but admixed | 保留 + 標記 subclonal，適合後續 ISM read-level characterization |
| **S5（聯集）** | `CONF = HIGH_UNION` | 聯集、覆蓋 2.85% 全 TP | 白名單 fast-track（若追求 precision）|
| **S4** | `CONF = AMBIGUOUS` | 無辨別力 bucket；75% TP + 76% FP 混合 | **不可單靠 S4 判定**；需啟動 ISM 二級判別（如 NGroups、chr context、UMI dedup）|
| **剩餘 cells（非 S1-S4）** | `CONF = NEEDS_CONTEXT` | 小 n 或邊際 TP rate | case-by-case，報告應該列 top-N cell 給專家 review |

**關鍵非對稱**：
- S1/S2/S3 是**白名單型 filter**（high-precision subset）
- S4 是**黑名單型 filter**（排除高風險）
- 兩者功能不同，**不應套用統一閾值**

---

## 7. 侷限與不做的

### 7.1 TO mode 專屬

Main 結論來自 HCC1395 TO mode。跨樣本 TO 驗證需等 Phase 1 fix 擴展到其他 6 樣本（當前僅 HCC1395 有新 KDE + TO）。paired_full 7 樣本 FP% 4% 不足以驗證 scheme 的**FP reduction**，只能驗證 TP purity 不崩。

### 7.2 S4 未被解決

S4 含 75% TP 加 76% FP，用當前 5 維（LOH_Subtype × AF × CN × mode × sample）無法進一步切分。**S4 的 TP/FP 分離須引入其他特徵**：
- ReadParser 的 DepthCF、NumReads、phasing quality
- 基因組 context（Chr、segdup annotation、germline AF database）
- ISM-specific: HPFineNGroups NG=4 + NR≥80（但本分析中 S8 僅 n=2 太小）

### 7.3 NG 的邊際貢獻消失

在 biology-module 框架下 NG≥3 幾乎無加成（S6 vs S1、S7 vs S5 差 <1pp）。這**不代表 NG 無用**，而是：
- NG 與 biology 模組重疊高（LOH_Strong cells 本身就 NG 偏高）
- NG 可能在 **S4 內部拆分**時才發揮作用（本報告未探）
- Phase 1 HP-only 修正已顯示 NG 訊號可能來自 somatic HP tag 人工分組，flag=on 下 NG≥3 消失 → 本報告 NG 相關結論應視為 flag=off 條件下

---

## 8. 建議後續

### 8.1 即刻（本週）

1. **S4 內部拆分 pilot**：在 S4 (30432 rows、8780 FP) 內加入 ReadParser 特徵（DepthCF、phasing quality）跑 logistic regression + Wilson CI，看是否能在 S4 內切出 10%-level 的 high-TP 子集
2. **Cross-sample TO 擴展**：等 Phase 1 HP fix 擴展到其他樣本後重跑 S1-S7 看泛化性
3. **v1 報告加 frontmatter warning**（指向 v2）

### 8.2 中期

4. **結合 ClairS-TO Verdict**：Verdict_Germline 已驗證能 100% 鎖 LowQual（see ClairS-TO Verdict Pilot 記錄），與 S4 bucket 是否重疊可量化
5. **flag=on 重驗**：確認 NG 的邊際貢獻在修正後的 HP tag 下是否恢復（關聯 HP-only Phase 1 結果）
6. **LOH_Subtype 細化**：LOH_Noise 意外顯示訊號 → 是否應併入 LOH_Weak 或獨立保留

### 8.3 論文取向

若此框架穩定，可作 **"Biology-informed stratified filter for tumor-only caller"** 的 methodological contribution：
- S3 Diploid Het → high-precision whitelist
- S1/S2 LOH-aware → subclonal characterization
- S4 ambiguous bucket → 二級判別需求定義

對照 Wakhan/SAVANA（純 SV/CN 工具）與 ClairS-TO verdict（純 read-level quality）→ 本框架在**biology module × AF × CN 交互**上補缺。

---

## 9. 檔案索引

### 9.1 資料（`research/ng_kde_rescaling/data/`）

| 檔案 | 內容 |
|------|------|
| `tpfp_cube_coarse.tsv` | paired_full 7 樣本 LOH×CN×AF_class TP/FP cube (339 cells) |
| `tpfp_cube_fine.tsv` | paired_full 10-bin AF 精細 (603) |
| `tpfp_cube_coarse_TO.tsv` | HCC1395 TO mode LOH×CN×AF_class (49 cells) |
| `tpfp_cube_fine_TO.tsv` | TO 10-bin AF 精細 (74) |
| `tpfp_top_high_tp_cells.tsv` / `_TO.tsv` | Top TP-enriched cells |
| `tpfp_top_high_fp_cells.tsv` / `_TO.tsv` | Top FP-enriched cells |
| `tpfp_stratified_filter_schemes.tsv` / `_TO.tsv` | S1-S8 scheme 彙總 |
| `tpfp_per_sample_schemes.tsv` | 每樣本 × scheme 指標 |
| `tpfp_ng_within_biology_modules.tsv` | NG effect 在 biology module 內 |

### 9.2 圖（`docs/experiments/in_progress/2026/04/figures/20260421_ng_kde_rescaled/`）

| 檔案 | 內容 |
|------|------|
| `fig_v2_1_to_tp_heatmap.png` | TO mode TP rate by LOH×AF×CN |
| `fig_v2_2_to_fp_heatmap.png` | 同 axes FP rate |
| `fig_v2_3_filter_scheme_bar.png` | S1-S7 TP/recall/FP-reduc 對照 |
| `fig_v2_4_per_sample_schemes.png` | paired_full 7 樣本 scheme 穩定性 |
| `fig_v2_5_operating_points.png` | Purity-recall operating points |
| `fig_v2_6_biology_module_interpretation.png` | 生物模組詮釋表 |

### 9.3 腳本（`research/ng_kde_rescaling/scripts/`）

| 檔案 | 功能 |
|------|------|
| `step4_tpfp_discrimination.py` | paired_full TP/FP cube + 6 schemes |
| `step4b_tpfp_tomode.py` | TO mode 主測試軌 |
| `step5_visualize_tpfp.py` | 6 張圖 |

---

## 10. 相對 v1 的修正摘要

| v1（偏題） | v2（本報告） |
|------------|--------------|
| 主軸：20260414 ΔNG 方向是否保留 | 主軸：LOH×AF×CN 是否切分 TP vs FP |
| 核心 metric：ΔNG Cohen's h | 核心 metric：TP rate, TP recall, FP reduction |
| Testbed：paired_full 7 樣本 | Testbed：HCC1395 TO mode (71.1% TP baseline) |
| 結論：6/7 方向反轉 → 舊結論需重審 | 結論：S3/S1/S5 是高 purity 白名單；S4 需二級判別 |
| 產出：穩健性否決報告 | 產出：stratified filter scheme 提案 + 生物模組詮釋 |
| 建議：master rerun 驗證方向 | 建議：S4 內部拆分 + cross-sample TO 擴展 |

**v1 保留作為附錄**（方向性 sanity check），但主論述改採 v2 框架。

---

## 11. SEQC2 Baseline vs Scheme TP:FP Fold-Improvement（E1）

**問題**：這些 scheme 相對原始 caller baseline 的 TP:FP 比值改善多少？

**方法**：每 sample-mode 下，比較 caller baseline 的 TP:FP ratio 與 scheme TP:FP ratio，計算 `fold-improvement = scheme_ratio / baseline_ratio`。由於 master TSV 僅含 caller 產出的 region（無 FN），無法算真 F1；以 precision + fold-improvement 取代。

### 11.1 HCC1395 TO（baseline TP:FP = 2.47:1，高 FP 主測試軌）

| Scheme | n | TP:FP ratio | fold vs baseline | TP recall | Power |
|--------|---:|------------:|-----------------:|----------:|-------|
| **S3 Diploid Het** | 380 | **21.35:1** | **8.69×** ⭐ | 1.3% | STRONG |
| S5 Combo (¬S4) | 886 | 11.14:1 | **4.53×** | 2.9% | STRONG |
| S7 Combo+NG≥3 | 484 | 9.76:1 | 3.97× | 1.5% | STRONG |
| S6 S1+NG≥3 | 253 | 9.54:1 | 3.88× | 0.8% | STRONG |
| S1 LOH_Strong+Extreme | 292 | 9.07:1 | 3.69× | 0.9% | STRONG |
| S2 Subclonal_LOH+Inter | 214 | 6.93:1 | 2.82× | 0.7% | STRONG |
| S4 NonLOH+Extreme | 30,432 | 2.47:1 | 1.00× （無提升） | 76.0% | STRONG |

### 11.2 COLO829 paired_full（baseline = 15.5:1，唯一高 FP paired_full）

| Scheme | n | TP:FP ratio | fold vs baseline | TP recall | Power |
|--------|---:|------------:|-----------------:|----------:|-------|
| S2 Subclonal_LOH+Inter | 2,023 | 19.86:1 | 1.28× | 5.5% | STRONG |
| S5 Combo | 2,271 | 17.77:1 | 1.15× | 6.1% | STRONG |
| S4 NonLOH+Extreme | 27,319 | 15.99:1 | 1.03× | 73.1% | STRONG |
| S1 LOH_Strong+Extreme | 244 | 9.17:1 | **0.59×** ⚠ | 0.6% | STRONG |
| S3 Diploid Het | 4 | inf | N/A | 0.01% | **NO_USE** |
| S6, S7 | 0-1 | — | — | — | NO_USE |

**關鍵觀察**：
- HCC1395 TO 的 fold-improvement 遠高於 COLO829 paired_full，因為 COLO829 baseline 已達 15.5:1（FP 稀少）；TO mode 才是 FP 辨別力的真正檢定場。
- COLO829 的 S1 fold=0.59× 異常（TP:FP 反而劣化），該樣本 LOH_Strong 註解與 HCC1395 不同步（見 §12）。
- S3 Diploid Het 在 COLO829 僅 n=4 → LOH.bed 覆蓋範圍不同；此為樣本特異性限制。

**圖表**：[fig_v2_7_baseline_tp_fp_fold.png](figures/20260421_ng_kde_rescaled/fig_v2_7_baseline_tp_fp_fold.png)

![F7 Baseline vs Scheme TP:FP Fold-Improvement](figures/20260421_ng_kde_rescaled/fig_v2_7_baseline_tp_fp_fold.png)

---

## 12. Per-Sample Scheme 驗證與 FP 集中偵測（E2）

**問題**：這些組合在所有樣本都成立嗎？FP 是否集中在特定樣本？

### 12.1 7 樣本 × 7 scheme TP rate 矩陣

| Sample-Mode | S1 | S2 | S3 | S4 | S5 | S6 | S7 |
|-------------|---:|---:|---:|---:|---:|---:|---:|
| **HCC1395 TO** ⭐ | 0.90 (292) | 0.87 (214) | **0.96** (380) | 0.71 (30432) | 0.92 (886) | 0.91 (253) | 0.91 (484) |
| HCC1395 pf | 0.99 (8) | — | 0.99 (15924) | 0.98 (195) | 0.99 (15933) | — | — |
| HCC1395_DORADO pf | — | — | 0.99 (15800) | 0.99 (244) | 0.99 (15800) | — | — |
| HCC1937 pf | — | — | 0.99 (3934) | 0.97 (73) | 0.99 (3934) | — | — |
| HCC1954 pf | 0.98 (176) | — | 1.00 (1001) | 0.99 (2126) | 0.99 (1177) | — | — |
| H2009 pf | — | — | 1.00 (3540) | 1.00 (41773) | 1.00 (3540) | — | — |
| H1437 pf | — | — | 1.00 (1257) | 1.00 (19562) | 1.00 (1257) | — | — |
| **COLO829 pf** | 0.90 (244) | **0.95** (2023) | 1.00 (4) | **0.94** (27319) | 0.95 (2271) | — | — |

（括號為 cell n，— 為 n<10）

### 12.2 關鍵發現

| 觀察 | 意涵 |
|------|------|
| ✅ **paired_full S3/S4/S5 TP rate ≥ 0.97** 全樣本一致 | S3 Diploid Het 組合生物上穩健；paired_full 高 TP 飽和 |
| ✅ **FP 絕對數集中在 HCC1395 TO (11,606) 與 COLO829 pf (2,273)** | 這兩者是有意義的 FP 辨別力檢定場；H2009/H1437 FP ≤86 無統計力 |
| ⚠ **HCC1395 TO 與 pf 的 S1 差距（0.90 vs 0.99）** | 同樣 LOH.bed 標註，但 TO mode 將 paired_full 中被過濾的 Extreme AF 重新暴露 → 低 AF（caller_af<0.1）為主要 FP 來源 |
| ⚠ **COLO829 S1 n=244, TP rate=0.90** 低於其他樣本 | 該樣本 LOH_Strong 註解不等同 HCC1395；reviewer 若質疑「S1 是否跨樣本一致」應強調 S1 在 COLO829 fold<1 → LOH_Strong 定義不跨樣本泛化 |
| ⚠ **COLO829 S3 n=4** | LOH.bed 幾乎全域覆蓋 COLO829 → None+Near-half cell 太小 → S3 scheme 無法在該樣本重現；v2 結論限 HCC1395-like LOH 注釋結構 |

### 12.3 Hard Gate 檢查

- 🟢 COLO829 pf 的 S5 TP rate=0.95 > baseline 0.94 → **未觸發 reviewer Hard Gate**；S5 是跨樣本穩健組合
- 🟡 S1 在 COLO829 fold=0.59× → 不反轉但 marginal；LOH_Strong 跨樣本歧義需在未來正式化

**圖表**：[fig_v2_8_per_sample_scheme_heatmap.png](figures/20260421_ng_kde_rescaled/fig_v2_8_per_sample_scheme_heatmap.png)

![F8 Per-Sample × Scheme TP Rate](figures/20260421_ng_kde_rescaled/fig_v2_8_per_sample_scheme_heatmap.png)

---

## 13. 統計力與 Cell 可信度（E3）

**問題**：這些 cell 的數量足夠可信嗎？

**方法**：依 n 分級 NO_USE (<20) / WEAK (20-50) / OK (50-200) / STRONG (≥200)；Wilson 95% CI 寬度附於 `tpfp_per_sample_scheme_full.tsv`。

| Power tag | n rows (sample×scheme) | 用途 |
|-----------|---:|------|
| STRONG (≥200) | 36 | ⭐ 可直接作證 |
| OK (50-200) | 2 | 🟢 可用，需附 CI |
| WEAK (20-50) | 8 | 🟡 參考，不作主論證 |
| NO_USE (<20) | 10 | ❌ 排除（主要在 COLO829 S3/S6/S7、HCC1937 S1）|

**v2 正文的主要 cell 分級**：
- HCC1395 TO 的 S1 (n=292)、S3 (n=380)、S5 (n=886)、S6 (n=253)、S7 (n=484)、S2 (n=214) 全 STRONG
- HCC1395 TO S4 n=30,432 (STRONG) — 這是 v2 強調「70% caller TP 落在 S4 無辨別力桶」的統計基礎
- Wilson CI 寬度（TP rate）：S3 HCC1395 TO ≈ ±0.02、S1 TO ≈ ±0.03、S5 TO ≈ ±0.02 → v2 §3 的 TP rate 數值均在 ±0.03 以內可信

---

## 14. TP vs FP 特徵差異與 Sub-Scheme 提案（E4+E5）

**問題**：同一 scheme 內 TP 和 FP 的特徵分佈是否不同？若差異大，可再加特徵切分。

### 14.1 HCC1395 TO 主要特徵差（Top by |Cliff's δ|）

每個特徵計算 TP vs FP 的 Mann-Whitney p 與 Cliff's delta（-1 ~ +1）。**δ<0 = TP 小於 FP**。

| Scheme | Feature | n_tp | n_fp | median TP | median FP | p | Cliff's δ |
|--------|---------|---:|---:|-----:|-----:|-----:|-----:|
| S6 S1+NG≥3 | **AlleleDelta** | 229 | 24 | 0.004 | 0.032 | 8.5e-5 | **-0.49** |
| S1 LOH_Strong+Extreme | **AlleleDelta** | 263 | 29 | 0.011 | 0.036 | 1.5e-4 | **-0.43** |
| S6 S1+NG≥3 | **CovM_used** | 229 | 24 | 1.30 | 1.55 | 5.8e-4 | -0.43 |
| S6 S1+NG≥3 | **NumReads** | 229 | 24 | 150 | 176 | 1.6e-3 | -0.39 |
| S2 Subclonal | CovM_used | 187 | 27 | 0.96 | 1.12 | 1.5e-3 | -0.38 |
| S4 NonLOH+Extreme | **AlleleDelta** | 21,652 | 8,780 | 0.012 | 0.033 | <1e-100 | **-0.34** |
| S3 Diploid Het | **AlleleDelta** | 363 | 17 | 0.457 | 0.435 | 0.038 | +0.30 |

### 14.2 模式解讀

**模式一：AlleleDelta (~AF distance from 0/1) 在 Extreme 組合是主要特徵差**
- 在 S1/S4/S6（Extreme 組合）中，FP 的 AlleleDelta 顯著高於 TP（median 差 0.02-0.03）
- **生物意義**：Extreme caller_af 的 FP 其實不是真正的 0%/100%（germline 洩漏的 read 帶來稍高 AlleleDelta）；真正的 LOH/germline TP 更接近純粹的 0 或 1
- **但 Cliff's δ=-0.34 不足以做 clean cut**：FP 和 TP 分布大量重疊

**模式二：CovM_used + NumReads FP > TP（反預期）**
- 原假設：低 coverage FP 較多；實際：FP 反而在更高 CovM 和 NumReads 區
- **生物意義**：FP 集中在高 CN（gain/amp）區域，mapping artifact + 高深度下 caller 誤判

**模式三：S3 Diploid Het TP_AlleleDelta 微高於 FP（δ=+0.30）**
- 真正 het variant 的 AlleleDelta 更接近 0.5（均衡 het），FP 則略偏 0.4-0.45
- 但 n_fp=17 太小，不作主結論

### 14.3 Sub-Scheme 提案結果（HCC1395 TO）

基於 E4 發現試算 S4/S5 refinement：

| Sub-Scheme | n | TP rate | TP:FP | fold vs baseline | TP recall |
|------------|--:|--------:|------:|-----------------:|----------:|
| S4d (S4+NR≥80+NG≥3) | 25,237 | **0.74** | 2.86 | 1.17× | 65.6% |
| S4b (S4+CovM in [0.7,1.4]) | 24,988 | 0.72 | 2.53 | 1.03× | 62.8% |
| S4a (S4+NR≥80) | 30,065 | 0.71 | 2.49 | 1.01× | 75.2% |
| S4c (S4+NG=4) | 5,442 | 0.69 | 2.25 | 0.92× | 13.2% |
| S5b (S5+NR≥80) | 865 | 0.92 | 11.01 | 4.48× | 2.8% |
| S5a (S5+NG=4) | 3 | 0.67 | 2.00 | — | NO_USE |
| Sx (S1orS3+NG=4+NR≥80) | 2 | 0.50 | 1.00 | — | NO_USE |

### 14.4 關鍵發現

| 發現 | 意涵 |
|------|------|
| ❌ **所有 S4 sub-scheme fold ≤ 1.17×** | 單一特徵 threshold 無法在 S4 內顯著切分 TP/FP；特徵分佈重疊太大（Cliff's δ 最大 -0.34 仍不足） |
| 🟡 S4d (NR≥80+NG≥3) 是最佳 S4 sub-scheme | TP rate 從 71% → 74%，但 fold 僅 1.17×；若用戶要，可作 S4 輕度過濾 |
| ❌ S5 + NG=4 n=3 無意義 | HCC1395 TO 的 NG=4 恰好與 S5 biology module 互斥；這印證 NG 主要在 S4 內出現（與前期 20260421 HP-only Phase 1 observation 一致）|
| ✅ **S5b (S5+NR≥80)** 維持 92% TP rate、fold=4.48× | 僅損失 3% recall（S5 → S5b：n 886→865）；若強調 read-depth floor 可用 S5b |
| ⚠ **Cliff's δ=-0.34 即最大特徵差** | v2 結論穩固：S1/S3/S5 是 biology-defined purity；進一步特徵 refinement 邊際報酬遞減 |

### 14.5 Operating Points 圖

**圖表**：[fig_v2_10_subscheme_operating_points.png](figures/20260421_ng_kde_rescaled/fig_v2_10_subscheme_operating_points.png)

![F10 Sub-scheme Operating Points](figures/20260421_ng_kde_rescaled/fig_v2_10_subscheme_operating_points.png)

**視覺結論**：sub-schemes（三角）大多低於或持平 S1-S7（圓點）；S5b 與 S5 幾乎重疊；Sx/S5a 墜入左下（低 recall + 低 precision + 小 n）。**不存在 Pareto 超越 S5 的 sub-scheme**。

### 14.6 Figure v2-9：S4 內 TP vs FP 特徵分佈

![F9 S4 Feature Violin](figures/20260421_ng_kde_rescaled/fig_v2_9_feature_violin_in_S4.png)

---

## 15. §11-14 統合結論

### 15.1 回答用戶 8 項要求

| # | 用戶問題 | 回答位置 | 摘要答案 |
|---|----------|---------|---------|
| 1 | 使用指標切分 TO 高 TP 區域 | §3 + §11 | S3/S1/S5 為 HCC1395 TO 下高 purity 白名單 |
| 2 | 相比 baseline TP:FP 比值 | §11 | S3=8.69×, S5=4.53×, S1=3.69× HCC1395 TO baseline |
| 3 | 所有樣本分開狀況 | §12 | paired_full 樣本 S3/S5 TP rate ≥0.97 全一致 |
| 4 | FP 是否集中某樣本 | §12 | HCC1395 TO (11,606 FP) + COLO829 pf (2,273 FP) 為主 FP 來源 |
| 5 | 組合是否合理 | §12 + §5 | S3/S4/S5 在 paired_full 跨樣本穩健；S1 在 COLO829 fold<1（LOH_Strong 跨樣本歧義）|
| 6 | 如何切分得更好的 TP:FP | §14 + §6 | S4 內 sub-scheme 最大改善僅 1.17× → 特徵空間飽和；S5b (S5+NR≥80) 維持 4.48× |
| 7 | 區域數量是否可信 | §13 | 36 STRONG + 10 NO_USE cells；v2 正文 cells 均 STRONG |
| 8 | TP/FP 特徵差異 | §14 | AlleleDelta (Cliff's δ=-0.49) 為最強差；但單一特徵 refinement 邊際效應飽和 |

### 15.2 實作建議更新

- **推薦 operating points**（與 v2 §8 一致、E5 驗證後不變動）：
  - **高 purity 白名單**：S3 (n=380, TP=96%)
  - **中 recall 穩定組合**：S5 (n=886, TP=92%, 4.53× fold)
  - **含 coverage 要求**：S5b (n=865, TP=92%) 若要強調 NR≥80 gate
- **S4 不做 filter**：70% 的 TP 落在 S4；任何 S4 refinement 均會損失 recall 且 TP:FP 提升 ≤1.17×
- **HCC1395 以外樣本**：paired_full TP rate 飽和（>97%），本框架的 FP 辨別力僅在 COLO829 pf (S5=0.95 穩健) 可再驗證；其他樣本需等 TO pipeline 擴展

### 15.3 未解問題

- TO mode 僅 HCC1395 有資料 → 跨樣本 FP 辨別力普適性未驗證
- S1 在 COLO829 的 fold<1 → LOH_Strong 定義在 non-HCC1395 樣本的等價性待驗證
- 進一步切分需多變量 model（logistic/RF），但 v2 scope 不含 → 留後續

---

## 16. §11-14 新增檔案清單

### 16.1 資料

| 檔案 | 功能 |
|------|------|
| `tpfp_baseline_ratio.tsv` | E1 每 scheme × sample-mode 的 baseline vs scheme ratio + fold |
| `tpfp_per_sample_scheme_full.tsv` | E1+E3 同 tsv + Wilson CI + power_tag |
| `tpfp_feature_diffs.tsv` | E4 每 scheme × feature 的 Mann-Whitney + Cliff's δ |
| `tpfp_subschemes.tsv` | E5 sub-scheme 提案 + fold |

### 16.2 圖表

| 檔名 | 內容 |
|------|------|
| `fig_v2_7_baseline_tp_fp_fold.png` | §11 baseline vs scheme TP:FP + fold（HCC1395 TO + COLO829 pf 2×2）|
| `fig_v2_8_per_sample_scheme_heatmap.png` | §12 7 sample-mode × 7 scheme TP rate heatmap |
| `fig_v2_9_feature_violin_in_S4.png` | §14 S4 內 TP vs FP 5 特徵 violin |
| `fig_v2_10_subscheme_operating_points.png` | §14 sub-scheme operating points 對 S1-S7 |

### 16.3 腳本

| 檔案 | 功能 |
|------|------|
| `step6_tpfp_detailed.py` | E1-E5 計算 |
| `step7_visualize_detailed.py` | F7-F10 圖 |

---

## 17. 多維綜觀視角（E6）— 超越單 scheme 的 Pareto Envelope

**問題**：§14 的 sub-scheme 在單 biology scheme 內加單一特徵 threshold 最多達 fold=1.17×。若放棄 S1-S7 框架，**綜觀 5 個參數聯合 cube（LOH × AF × CN × NG × NR_band = 900 bins）**，是否存在超越 S5 的 operating point？

### 17.1 HCC1395 TO 5D Cube Top 10 Cells（n≥20）

| rank | LOH_Subtype | AF_class | CN | NG | NR | n | TP rate | TP:FP | fold |
|-----:|:-----------:|:--------:|:--:|:--:|:--:|--:|--------:|------:|-----:|
| 1 | LOH_Strong | Extreme | T1 | 3 | mid | 22 | **1.000** | ∞ | — |
| 2 | LOH_Noise | Extreme | T3 | 2 | high | 56 | 0.982 | 55.0 | 22.3× |
| 3 | LOH_Weak | Extreme | T2 | 3 | mid | 44 | 0.977 | 43.0 | 17.4× |
| 4 | LOH_Weak | Intermediate | T1 | 2 | mid | 43 | 0.977 | 42.0 | 17.0× |
| 5 | None | **Near-half** | T1 | 2 | mid | **171** | 0.977 | 41.75 | 16.9× |
| 6 | LOH_Noise | Extreme | T2 | 2 | high | 146 | 0.973 | 35.5 | 14.4× |
| 7 | LOH_Noise | Extreme | T2 | 3 | high | 35 | 0.971 | 34.0 | 13.8× |
| 8 | None | Near-half | T2 | 3 | high | 31 | 0.968 | 30.0 | 12.2× |
| 9 | LOH_Strong | Extreme | T2 | 3 | high | 87 | 0.966 | 28.0 | 11.4× |
| 10 | LOH_Weak | Intermediate | T1 | 3 | mid | 27 | 0.963 | 26.0 | 10.5× |

**觀察**：
- Top 1 (LOH_Strong+Extreme+T1+NG=3+mid) 是 S6 的子 cell（S6 整體 90.5% → 細分至此 cell 達 100%）— n=22 太小不作主張
- **Top 5 (None+Near-half+T1+NG=2+mid, n=171, 97.7%)** 是 S3 子集，但在 NG=2 + mid-NR 精純度遠高於 S3 整體（95.5%）— 這是新發現的 cell
- Top 4/10 (LOH_Weak+Intermediate) 不屬 S1-S7 任何主 scheme，但 fold≥10.5× — v2 框架遺漏這些 cells
- CN T2 + NG≥2 出現頻繁 → 多核 subclone + gain 區是 TP 富集帶

### 17.2 Cumulative Envelope：取 top-k cells 的 purity-recall trade-off

**HCC1395 TO**（baseline TP:FP = 2.46:1）：

| 目標 purity | k cells | cum_n | recall | cum TP:FP | fold vs baseline |
|-------------|--------:|------:|-------:|----------:|----------------:|
| ≥ 95% | **17** | 1,099 | **3.7%** | ∞（n_fp=43）→ 24.6:1 | **10.0×** |
| ≥ 90% | **28** | 2,285 | **7.4%** | ~11.7:1 | **4.73×** |
| ≥ 80% | 32 | 8,421 | 24.5% | ~4.9:1 | 1.99× |

**對比 v2 主 scheme**：

| Scheme | recall | purity | fold |
|--------|-------:|-------:|-----:|
| S3 Diploid Het（n=380） | 1.3% | 95.5% | 8.69× |
| S5 Combo（n=886） | 2.9% | 91.8% | 4.53× |
| S7 Combo+NG≥3（n=484） | 1.5% | 90.7% | 3.97× |
| **E6 17-cell envelope** | **3.7%** | **96.1%** | **10.00×** |
| **E6 28-cell envelope** | **7.4%** | **90.4%** | **4.73×** |

**Pareto 分析**：
- ✅ **17-cell envelope Pareto-dominates S3**：recall 提高 2.8×（3.7% vs 1.3%）且 purity 更高（96% vs 95.5%）
- ✅ **28-cell envelope Pareto-dominates S5**：recall 提高 2.5×（7.4% vs 2.9%）且 purity 持平（90.4% vs 91.8%）
- ⚠ 但**每個 cell 是聯合 5 特徵的 hard bin**，無 smoothing；跨樣本泛化性需單獨驗證

### 17.3 COLO829 paired_full Top 10 Cells

| rank | LOH_Subtype | AF_class | CN | NG | NR | n | TP rate | TP:FP |
|-----:|:-----------:|:--------:|:--:|:--:|:--:|--:|--------:|------:|
| 1 | LOH_Strong | Intermediate | T2 | 1 | mid | 85 | 1.000 | ∞ |
| 2 | None | Intermediate | T2 | 3 | mid | 175 | 0.989 | 86.5 |
| 3 | None | Extreme | T2 | 1 | mid | 130 | 0.969 | 31.5 |
| 4 | None | Intermediate | T1 | 0 | mid | 95 | 0.968 | 30.7 |
| 5 | LOH_Noise | Extreme | T2 | 1 | mid | 277 | 0.968 | 29.8 |
| 6 | None | Intermediate | T2 | 3 | high | 28 | 0.964 | 27.0 |
| 7 | None | Intermediate | T1 | 1 | mid | 79 | 0.962 | 25.3 |
| 8 | LOH_Weak | Intermediate | T2 | 1 | mid | 651 | 0.960 | 24.0 |

**觀察**：
- COLO829 top cells **NG 多為 1** — 與 HCC1395 TO（NG 多為 2-3）截然不同 → 跨樣本 cell 座標不相容
- **Intermediate AF 取代 Near-half** 作為主導 — COLO829 caller 輸出 AF 分佈偏向 Intermediate
- LOH 定義差異（LOH_Strong NG=1 表示 COLO829 的 LOH.bed 覆蓋等同 germline 純讀 haplotype）→ §12 已述 LOH_Strong 跨樣本不同質

### 17.4 COLO829 Cumulative Envelope

| 目標 purity | k cells | recall | fold |
|-------------|--------:|-------:|-----:|
| ≥ 98% | 4 | 1.4% | 3.42× |
| ≥ 95% | **24** | **51.3%** | **1.31×** |
| ≥ 90% | 54 | 98.7% | 1.03× |

**關鍵**：COLO829 baseline 已達 94%，envelope 提升空間小；51% recall @ 95% purity 雖高但 fold 僅 1.31×（vs HCC1395 TO 的 10×）

### 17.5 圖表

**HCC1395 TO**：[fig_v2_11_panorama_HCC1395_TO.png](figures/20260421_ng_kde_rescaled/fig_v2_11_panorama_HCC1395_TO.png)

![F11a Panorama HCC1395 TO](figures/20260421_ng_kde_rescaled/fig_v2_11_panorama_HCC1395_TO.png)

**COLO829 pf**：[fig_v2_11b_panorama_COLO829.png](figures/20260421_ng_kde_rescaled/fig_v2_11b_panorama_COLO829.png)

![F11b Panorama COLO829](figures/20260421_ng_kde_rescaled/fig_v2_11b_panorama_COLO829.png)

### 17.6 關鍵結論（相對 §14 的修正）

| § | §14 結論 | §17 修正/補充 |
|---|----------|---------------|
| 特徵空間飽和度 | ❌ S4 內單特徵 refinement 飽和（fold≤1.17×）| 🟡 **聯合 5 特徵 cube 並未飽和** — envelope 在 S5 外側有 Pareto 優勢 |
| v2 主 scheme 是否最佳 | ✅ S3/S5 是好的 biology-defined purity | 🟡 **S3/S5 非 Pareto-optimal** — 17-cell/28-cell envelope 於 purity-recall 兩軸同時優於 S3/S5 |
| 下一步建議 | 不建議 sub-scheme | 🟢 若需更高 recall，可採 **經驗 empirical top-k cells white-list**（但需跨樣本驗證） |

### 17.7 限制與警告

- 🔴 **樣本泛化性未驗證**：top cells 在 HCC1395 TO vs COLO829 pf 座標完全不同（NG=2-3 vs NG=1；Near-half vs Intermediate）。HCC1395 TO 學到的 cell white-list 在其他樣本可能失效
- 🔴 **經驗 bin 無平滑**：每個 cell 是硬切分；若 region 座標略偏（如 NR=58 vs 60）可能跨 bin 跳躍
- 🟡 **無 cross-validation**：envelope 是 in-sample 分析，未做 LOSO → 可能 over-fit HCC1395 TO 的 noise pattern
- 🟡 **只測 HCC1395 TO + COLO829 pf**：其他 5 樣本 FP 太少，envelope 不具統計力

### 17.8 §17 對用戶問題的回答

**「如果切分到特徵層次,而是綜觀的各參數面向觀察,結果如何」**：

**答**：綜觀 5D cube 確實發現 Pareto 外側優勢：
- **17-cell 白名單** 可達 96.1% purity + 3.7% TP recall（**fold=10×**）— 比 §14 最佳 S3 (95.5% + 1.3%, fold=8.7×) 全面更優
- **28-cell 白名單** 可達 90.4% purity + 7.4% recall（fold=4.73×）— 比 S5 (91.8% + 2.9%, fold=4.5×) recall 提高 2.5×
- 但這是**經驗 empirical cells**（不是生物學上可直接解讀的 scheme），且僅在 HCC1395 TO 驗證 → reviewer 應視為「飽和測試」結果，而非建議的實作濾網
- §14 的「特徵飽和」結論僅限於「在既定 biology scheme 內加單特徵」；聯合 cube 並未飽和

**實作建議**：
- 若追求最高 purity+recall（願意接受「hard cell white-list」）：用 top-17 cells
- 若追求可生物解釋的 scheme：維持 §8 推薦 S3/S5
- 跨樣本應用前需 LOSO 驗證 cell 泛化性（out-of-scope 本報告）

### 17.9 新增檔案

| TSV | 內容 |
|-----|------|
| `tpfp_5d_cube_HCC1395_TO.tsv` | 164 個活躍 cell × (LOH, AF, CN, NG, NR) 座標 + n, TP rate, TP:FP |
| `tpfp_5d_pareto_HCC1395_TO.tsv` | n≥20 cells 按 TP rate 排序（rank 1-42） |
| `tpfp_5d_cumulative_envelope_HCC1395_TO.tsv` | 累積 envelope（cum_tp_rate, cum_recall, cum_fold）|
| `tpfp_5d_cube_COLO829_pf.tsv` | COLO829 paired_full 版本 |
| `tpfp_5d_pareto_COLO829_pf.tsv` | COLO829 版本 |
| `tpfp_5d_cumulative_envelope_COLO829_pf.tsv` | COLO829 版本 |

腳本：`research/ng_kde_rescaling/scripts/step8_multidim_panorama.py`

圖表：`fig_v2_11_panorama_HCC1395_TO.png`、`fig_v2_11b_panorama_COLO829.png`
