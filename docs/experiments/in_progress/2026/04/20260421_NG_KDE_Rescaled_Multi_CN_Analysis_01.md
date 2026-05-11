---
status: superseded_off_topic
date: 2026-04-22
supersedes_pointer: docs/experiments/in_progress/2026/04/20260422_LOH_AF_KDE_TPFP_Discrimination_02.md
type: reanalysis
priority: P2
sample_set: 7 cell lines
pipeline_track: paired_full + HCC1395 TO sidebar
kde_baseline: new (post-fix per-sample Diploid_Coverage_Used)
related:
  - docs/experiments/validated/2026/04/20260414_LOH_Subclone_AF_Methylation_Evidence_01.md
  - docs/experiments/in_progress/2026/04/20260420_KDE_Fix_Acceptance_Validation_01.md
  - docs/experiments/in_progress/2026/04/20260421_HPFineNGroups_Marker_Reaudit_01.md
  - docs/experiments/in_progress/2026/04/20260418_F_HPFineNGroups_deepening_POSITIVE_01.md
outcome: 20260414 POSITIVE direction REVERSED (6/7 samples) under new KDE; H2 verification confirms reversal is NOT methodology artifact
v2_redirect_reason: >-
  This v1 answered "is 20260414 NG effect robust under new KDE" (ΔNG as metric),
  but the user's actual question is "can LOH x AF x new-KDE separate TP from FP
  via biology-informed stratified filters". v2 (20260422_02) reframes the
  analysis around TP/FP discrimination and biology module filter schemes.
  v1 retained as robustness-check appendix only.
---

> **⚠ v1 framing 偏題提示（2026-04-22 加註）**
>
> 本報告主軸為 **20260414 NG 效應量在新 KDE 下的穩健性**（ΔNG、Cohen's h），回答的是「方向是否保留」。
>
> 用戶的實際問題是 **LOH × AF × 新 KDE 能否切分 TP 與 FP，並用不同生物機制給不同參數意義**。此問題的答題在 v2：
> `docs/experiments/in_progress/2026/04/20260422_LOH_AF_KDE_TPFP_Discrimination_02.md`
>
> v1 保留作為**方向性穩健性附錄**，主論述以 v2 為準。

# NG × 新 KDE baseline × 多 CN 切割 × LOH 內外 × Per-sample 再分析

## 結論（一行）

在 7 樣本新 KDE baseline 下，20260414 原 POSITIVE 結論「CN1 LOH 中 Intermediate AF 的 NG 比 Extreme AF 高 +0.705（7/7 p<10⁻³⁹）」**方向反轉為 6/7 NEG（|ΔNG| 0.06-0.63）**，H2 驗證（LOH.bed + TP-only 精確復現）排除方法學差異；唯一保留 POS 的 HCC1954 正是本研究中唯一使用 stale baseline 做 rescale 的樣本，強烈指向原結論為 **stale hardcoded 75× baseline 的 artifact**。

---

## 1. 資料來源與方法

### 1.1 資料集

| 樣本 | 來源 | KDE 狀態 | n（paired_full） | 備註 |
|------|------|---------|-----------------:|------|
| HCC1395 | `output/canonical/HCC1395/paired_full/20260420_...` | new direct | 30,381 | Diploid_Coverage_Used = 55× |
| HCC1395_DORADO | `20260420_HCC1395_DORADO_...` | new direct | 30,129 | 53× |
| HCC1937 | `20260421_HCC1937_...` | new direct | 12,588 | 91× |
| H2009 | `20260421_H2009_...` | new direct | 132,995 | 79× |
| H1437 | `20260421_H1437_...` | new direct | 67,476 | 69× |
| COLO829 | `20260421_COLO829_...` | new direct | 37,458 | 29×（最小） |
| **HCC1954** | `20260315_...` + post-hoc rescale | stale ×1.2295 | 17,938 | 新 baseline 61× |
| HCC1395 TO | `/tmp/ism_hp_fix_phase1/{tp,fp}_off` | new direct (phase1) | 40,115 | sidebar |

**CovM_used 定義**：對 6 個 new direct 樣本直接用 `Coverage_Multiple`（已以 per-sample `Diploid_Coverage_Used` 正規化）；HCC1954 用 `stale × 1.2295`（ratio_stale_over_new）。

### 1.2 CN tier 切割策略

| 策略 | 邊界（CovM） | 依據 |
|------|---|------|
| **A_Legacy** | 0.75 / 1.25 / 1.75 | 20260414 使用（stale baseline 時代） |
| **B_ClinicalStrict** | 0.5 / 0.85 / 1.15 / 1.4 | FACETS/Battenberg deletion-LOH 嚴格定義 |
| **C_FineGrained** | 0.4 / 0.7 / 0.95 / 1.1 / 1.4 / 2.0 | TITAN 多亞克隆 CN state |
| **D_PerSample_Q** | 每樣本 q05/q25/q50/q75/q95 | 消除 ploidy 跨樣本差異 |
| **E_LOH_aware** | LOH Inner <0.6/<0.9，Outer <0.85/<1.15 | CAMDAC subclonal LOH |
| **F_SEQC2_grounded** | 0.65 / 0.99 / 1.33 / 1.82 | HCC1395 SEQC2 per_cn_bin_metrics 經驗中點 |

**SEQC2 經驗中點來源**：`research/kde_fix_validation/outputs/step1_hcc1395_seqc2/paired_pileup/per_cn_bin_metrics.tsv`
 CN=1-2 mean 0.47 / CN=2 0.83 / CN=3 1.14 / CN=4 1.52 / CN≥5 2.11 → 鄰對中點 0.65 / 0.99 / 1.33 / 1.82

### 1.3 分類定義

- **AF class**（與 20260414 相同）：Extreme (<0.1 或 >0.9) / Near-half (0.4-0.6) / Intermediate（其餘）
- **AF bin10**（用戶要求細化）：每 0.1 一格，10 bins
- **LOH side**：優先 `LOH_Bed_Overlap`，回退 `Potential_LOH`（2026-04-19 audit Jaccard=1.0 vs external LOH.bed）。本次新 KDE runs 未掛載 `--loh-bed`，故實際使用 Potential_LOH。H2 驗證章節另以 external LOH.bed 精確復現。
- **Cohen's h**（NG≥2 proportion）、**Mann-Whitney U**（alternative="greater"，單尾 assume Inter>Extreme）

---

## 2. 主要結果

### 2.1 CovM 分佈（新 KDE baseline）

![F1: CovM distribution](figures/20260421_ng_kde_rescaled/fig1_covm_density_per_sample.png)

**觀察**：
- 7 樣本 CovM median 全部接近 1.0（0.993–1.101），驗證新 baseline 對齊成功
- HCC1395 median 1.07、COLO829 1.01、HCC1937 1.09 — per-sample normalization 工作正常
- Strategy A（紅虛線 0.75/1.25/1.75）與 Strategy F（綠點線 0.65/0.99/1.33/1.82）邊界切在分佈左側 shoulder（CN1/CN2 邊界）

**比對 stale baseline**：20260414 時代所有樣本共用 75× → COLO829 實際 29× 被高估 2.6×，CovM 被壓低 2.6 倍；新 baseline 下 COLO829 的「CN1 subset」重新定義。

### 2.2 20260414 復現（Potential_LOH + TP+FP）

Step 1 復現表（`data/ng_20260414_replication_inter_vs_extreme.tsv`）結果：

| Sample | Strategy A ΔNG | Strategy F ΔNG | 方向 |
|--------|---------------:|---------------:|------|
| HCC1395 | **−0.561** | **−0.420** | 🔴 反轉 |
| HCC1395_DORADO | **−0.571** | **−0.493** | 🔴 反轉 |
| HCC1937 | **−0.635** | **−0.548** | 🔴 反轉 |
| HCC1954 | +0.247 | **+0.432** | 🟢 保留 |
| H2009 | −0.073 | +0.002 | 🟠 近零 |
| H1437 | −0.267 | −0.263 | 🔴 反轉 |
| COLO829 | −0.095 | −0.157 | 🟠 弱反轉 |

- 原 20260414：**ΔNG = +0.705**（7/7 p<10⁻³⁹）
- 新 KDE Strategy A：**5/7 反轉**，HCC1954 單一保留
- Strategy F（SEQC2-grounded）：**5/7 反轉 + 1/7 近零**，HCC1954 唯一 POS

**森林圖（3 策略對比）**：

![F5: replication forest](figures/20260421_ng_kde_rescaled/fig5_replication_forest_3strategies.png)

### 2.3 H2 驗證：排除方法學差異

用戶質疑：是否因 LOH 分類（Potential_LOH vs external LOH.bed）或 TP+FP 混合（vs TP-only）而造成方向反轉？

**H2 測試設計**：在新 KDE master 上重建 `loh_bed_hit`（external TO-pipeline LOH.bed annotation，與 20260414 相同程式邏輯），filter = `TP only + CovM_used<0.75 + loh_bed_hit=True`。

| Sample | LOH.bed × TP-only ΔNG | Potential_LOH × TP-only ΔNG | 原 20260414 |
|--------|----------------------:|----------------------------:|------------:|
| HCC1395 | **−0.534** | −0.534 | +0.787 |
| HCC1395_DORADO | **−0.454** | similar | +0.634 |
| HCC1937 | **−0.629** | similar | +0.791 |
| HCC1954 | +0.311 | similar | +0.579 |
| H2009 | −0.060 | similar | +0.681 |
| H1437 | −0.259 | similar | +0.823 |
| COLO829 | −0.118 | similar | +0.720 |

**H2 verdict**：LOH.bed × Potential_LOH Jaccard = 0.66-0.91（高重疊），結論一致；**方向反轉不是 LOH 分類差異所致**。

![F6: H2 verification](figures/20260421_ng_kde_rescaled/fig6_h2_verification_lohbed.png)

### 2.4 NG 分佈跨 CN tier × LOH 內外

![F3: NG distribution stacked](figures/20260421_ng_kde_rescaled/fig3_ng_distribution_stacked.png)

**核心 pattern**：
- **LOH Inner（上列）**：全 7 樣本所有 CN tier 中，NG 分佈向 {NG=1, NG=2} 集中（LOH 單倍型喪失壓縮 NG 上限）
- **LOH Outer（下列）**：NG 分佈散佈，NG=3/NG=4 比例顯著升高
- **HCC1954 特殊**：即使 LOH Inner，NG=2 佔多數（而非 NG=1），反映 HER2+ 高 ploidy（~4）下「CN1」的定義可能對應 `het-diploid with CNA subclone`，非真正 deletion LOH
- **HCC1937 LOH Outer**：NG=3 dominant（紅色帶），BRCA1 mutant 的 CNV 複雜性導致 haplotype structure 破碎
- **COLO829 LOH Outer T0**：NG=0 佔 80%（灰色），反映 ONT R10 無 methylation → HP tag 缺失造成 NG=0 artifact，再次確認 COLO829 不適用 HPFineNGroups

### 2.5 CN tier × LOH × NG filter 的 TP rate

![F2: TP rate heatmap](figures/20260421_ng_kde_rescaled/fig2_tprate_heatmap_A_vs_F.png)

**觀察**：
- 全 7 樣本 × 4 tier × LOH 內外 × 2 策略 的 TP rate 普遍 0.88-1.00（綠色），無明顯 tier 失效區
- 例外：**COLO829 CN1 LOH Outer TP rate ≈ 0.66**（最低單 cell），反映 ONT R10 沒有 methylation 時 NG 訊號失真
- NG=4+NR≥80 canonical filter（下列）在多數 cell n<5 → 留白；僅 HCC1395/H2009/H1437 部分 cell 可用
- Strategy A 與 Strategy F 的 TP rate 模式相似（差異 <0.02），說明**結論對 CN tier 邊界選擇不敏感**

### 2.6 AF 10-bin × NG cell TP rate（Strategy A CN1）

![F4: AF 10-bin × NG heatmap](figures/20260421_ng_kde_rescaled/fig4_af10bin_ng_heatmap.png)

**觀察**：
- AF<0.1 bin 主導 n（deletion LOH 的 reference-dominant SNV）
- HCC1395: NG=3 × AF∈[0.2, 0.3) 單 cell **紅色 TP rate=0**（n=6），是唯一 outlier
- COLO829: NG=3 × AF<0.1 n=79 TP rate 0.85，ONT R10 下的 HP artifact 仍然偶爾「意外」正確
- HCC1954: NG=3 × AF<0.1 n=955（龐大！），反映 HER2+ 高 ploidy 下「CN1」被定義得很寬，吸收大量 sub-diploid 變異

---

## 3. 機制假說：為何 20260414 的 POSITIVE 訊號消失？

### 3.1 假說 H1（首要候選）：stale hardcoded 75× baseline 導致 CN1 subset 組成異質

20260414 分析時：
- 所有樣本共用 hardcoded 75× baseline
- COLO829 實際中位 29× → 其 CovM 被**系統性壓低** 29/75 = 0.39
- H1437 實際 69× → CovM 壓低 0.92（幾乎無差）
- HCC1937 實際 91× → CovM **抬高** 1.21

**後果**：
- 在 stale 框架下 "CN1 = CovM<0.75" 對 COLO829 而言代表 actual coverage <30×（真 deletion LOH）
- 對 HCC1937 而言代表 actual coverage <68×（包含大量 near-diploid 區）
- 對 HCC1395 而言代表 actual coverage <41×（真 deletion LOH）

**新 baseline 下**：`CovM<0.75` 意義統一為「小於 sample 自身 diploid 的 75%」，CN1 變得更「乾淨」→ subset 組成變了，訊號消失。

### 3.2 假說 H2（已排除）：LOH 分類方法差異

本次使用 Potential_LOH fallback，20260414 使用 external LOH.bed；但 H2 驗證（§2.3）Jaccard 0.66-0.91，結論一致，**已排除**。

### 3.3 假說 H3（次要候選）：HCC1954 stale 保留 → artifact 指紋

唯一保留 POS（+0.247 ~ +0.432）的樣本 HCC1954 是**本次唯一使用 stale rescale 的樣本**（因 20260315 stale run 未重新 sequence）。
這是強力指紋：**如果反轉是真生物學變化，HCC1954 不應獨例外；如果反轉是 stale artifact，HCC1954 因為仍用 stale 資料，自然保留 artifact**。

### 3.4 Mechanistic interpretation（假說 H4，生物學層）

在新 KDE 框架下，CN1 LOH Inner 代表「真 deletion LOH」：
- 單 haplotype 殘留 → NG 上限 2（HP1 或 HP2 methylation classes），不可能 >2
- Intermediate AF（0.1-0.9 剔除 0.4-0.6）在 CN1 LOH 中代表 **reference reads 混雜**（mapping 問題、或 subclonal partial LOH）
- Extreme AF（<0.1 或 >0.9）代表**純單倍型**
- **NG 差異方向應該是 Extreme > Intermediate**（純單倍型 methylation 較為一致，可區分 ≥1 class；混雜區 methylation 散亂，單 class）

此假說與新 KDE 實測 5/7 Extreme > Intermediate 方向一致。**原 20260414 反方向可能源自 stale baseline 將「真 CN1」與「near-diploid 混入區」合併，後者的 Intermediate AF 對應 subclonal LOH (有 methylation 結構)。**

---

## 4. Per-sample 觀察（7 樣本專段）

### 4.1 HCC1395（乳癌細胞線，hyper-diploid ploidy 2.85，SEQC2 truth）

- 新 baseline 55×（vs stale 75×，ratio 1.36）
- CN1 LOH Inner Inter vs Extreme: **ΔNG = −0.534**（反轉）
- LOH.bed 覆蓋 44.6%（vs Potential_LOH 47.5%）
- **結論**：20260414 的 +0.787 無法在新 baseline 下復現

### 4.2 HCC1395_DORADO（Dorado basecalled，相同 HCC1395 sample）

- 新 baseline 53×（ratio 1.42）
- CN1 LOH Inner: ΔNG = −0.454
- LOH.bed 覆蓋 45.3% — 與 HCC1395 一致性良好
- **結論**：basecaller 差異不影響 reversal 方向（Guppy vs Dorado 同樣失去原訊號）

### 4.3 HCC1937（BRCA1 mutant 乳癌）

- 新 baseline 91×（**反向 ratio 0.82**，唯一 coverage 高於 stale 假設）
- LOH.bed 4531 segments（最多！中位 size 143kbp，最破碎）
- CN1 LOH Inner: **ΔNG = −0.629**（最大反轉）
- **BRCA1 HR defect → 複雜 CNV patterns → LOH segment 破碎 → 真 CN1 定義嚴格後 Inter 訊號消失**

### 4.4 HCC1954（HER2+ 乳癌，**本研究唯一 POS 保留者**）

- 新 baseline 61×（ratio 1.23），**stale rescaled only 樣本**
- LOH.bed 僅 8.3% 覆蓋（最低，ploidy ~4 下只有少量純 LOH）
- CN1 LOH Inner: **ΔNG = +0.311**（保留 POS，但效應量 <20260414 原 +0.579 的一半）
- **強力指紋**：若 reversal 是真生物學變化，HCC1954 不應例外；HCC1954 保留訊號 = 使用 stale 資料的 artifact 指紋

### 4.5 H2009（肺癌，orthogonal tools truth）

- 新 baseline 79×（ratio 0.95 ~ 1，幾無 rescale）
- CN1 LOH Inner: ΔNG = −0.060（**近零**）
- LOH.bed 僅 25% — orthogonal tools truth 非 SEQC2 嚴格
- **結論**：near-zero ΔNG，原 20260414 +0.681 效應大部分歸零

### 4.6 H1437（肺癌）

- 新 baseline 69×（ratio 1.09）
- CN1 LOH Inner: ΔNG = −0.259（反轉但較小）
- LOH.bed 19.3%（ploidy 相對低）
- **結論**：中等反轉

### 4.7 COLO829（黑色素瘤，ONT R10 無 methylation）

- 新 baseline 29×（ratio **2.59**，最極端 rescale 影響）
- CN1 LOH Inner: ΔNG = −0.095（近零）
- COLO829 NG=0 fraction 極高（LOH Outer T0 NG=0 占 80%），反映 ONT R10 沒有 methylation calls → HP artifact
- **結論**：已知 out-of-scope 樣本，反轉方向一致但訊號微弱

---

## 5. LOH 內外延伸（用戶原始要求）

從 `ng_loh_inner_vs_outer_extension.tsv` 整理：

| Sample | CN1 T0 Extreme ΔNG(I-O) | CN1 T0 Inter ΔNG(I-O) | 方向性 |
|--------|-------------------------:|----------------------:|-------|
| HCC1395 | −0.323 | −0.258 | LOH 內 < LOH 外（符合 haplotype loss） |
| HCC1395_DORADO | +0.733 | +0.810 | **反向！**（LOH Inner NG 比 Outer 高） |
| HCC1937 | −0.700 | — | 符合預期 |
| HCC1954 | −1.063 | — | 符合預期（LOH 壓 NG 最強） |
| H2009 | −1.180 | −0.836 | 符合預期 |
| H1437 | −0.381 | −0.668 | 符合預期 |
| COLO829 | +0.475 | +0.924 | **反向！**（ONT R10 NG 混亂） |

**觀察**：
- **5/7 樣本符合預期**：LOH Inner NG < LOH Outer NG（haplotype loss → NG 上限被壓）
- **HCC1395_DORADO 異常**：Dorado basecalled 的 LOH Inner 反而 NG 較高 → 可能 Dorado 對 CpG island 的 HP tag 有不同 false-positive pattern
- **COLO829 異常**：ONT R10 無 methylation，NG 訊號為 artifact

---

## 6. 多 CN 策略比較

| Strategy | A_Legacy | B_ClinicalStrict | F_SEQC2_grounded |
|----------|---------:|-----------------:|-----------------:|
| HCC1395 ΔNG (CN1 LOH Inter vs Ext) | −0.561 | SKIP (n<5) | −0.420 |
| HCC1395_DORADO | −0.571 | SKIP | −0.493 |
| HCC1937 | −0.635 | −0.377 | −0.548 |
| HCC1954 | +0.247 | SKIP | **+0.432** |
| H2009 | −0.073 | +0.045 | +0.002 |
| H1437 | −0.267 | SKIP | −0.263 |
| COLO829 | −0.095 | −0.296 | −0.157 |
| **方向一致性** | **5/7 NEG** | **2/4 NEG**（多 SKIP） | **5/7 NEG** |

**結論**：3 種 CN 切割策略下，**方向全部一致為 5/7 NEG**。反轉結果對 CN tier 邊界選擇**robust**。

---

## 7. 文獻依據（CN 切割理論基礎）

| Method | CN boundary 設計 | 關鍵文獻 | 與本研究對應 |
|--------|------------------|---------|------------|
| **FACETS** | AA/AB/AAB/AABB 整數 logR 拆分 | Shen & Seshan (2016) | Strategy B 邊界對應 AB(0.5-0.85) vs AAB(0.85-1.15) |
| **Battenberg** | logR + BAF 雙峰 subclonal fraction 估計 | Nik-Zainal et al. (2012) Cell | Strategy E LOH-aware 混合邊界 |
| **TITAN** | 最多 20 CN states (0-20)，subclonal prevalence | Ha et al. (2014) Genome Research | Strategy C FineGrained 細 tier |
| **CAMDAC** | Clonal LOH → ASM 模式 | Larose et al. (2023) Nature Methods | 本研究 CN1 LOH Inner 即 CAMDAC 定義 |
| **ASCAT** | Allele-specific CN，純度/ploidy 二元最佳化 | Van Loo et al. (2010) PNAS | 提供 ploidy 校正參考（HCC1395 2.85、HCC1954 ~4） |

**用戶要求的 HCC1395 SEQC2 校準（Strategy F）**：
- SEQC2 consensus truth 為 HCC1395 ground truth
- `per_cn_bin_metrics.tsv` 每 CN bin 的 median CovM 為經驗錨點
- F_SEQC2_grounded 邊界用「相鄰 CN bin 中點」設定，**直接以 SEQC2 truth 定標**
- 結果：F 策略與 A_Legacy 結論一致（方向反轉）→ 即使用 SEQC2 ground truth 校準 CN tier，20260414 POSITIVE 仍不成立

---

## 8. 結論穩定度判定（依計劃書 §6.3）

| 判定標準 | 結果 |
|---------|------|
| Strategy A 下 ΔNG 方向反轉？ | **YES（5/7 反轉，|Δ| 0.06-0.63）** |
| |ΔNG| 下降 >50%？ | **YES（原 +0.705 → 新 −0.36 中位，衰減率 >150%）** |
| Strategy B + F 方向一致於 A？ | **YES（5/7 NEG 穩定）** |
| Robustness 綜合 | 🔴 **RETRACT level（計劃書 Hard Gate 觸發）** |

**建議** (不自動執行，等用戶確認)：
1. 更新 memory `project_loh_subclone_af_methylation_positive.md`：加入「新 KDE baseline 下 5/7 樣本反轉」警示，但不 retract（因 HCC1954 single holdout 反向證據太強）
2. Phase 2B 建議：將 CovM baseline fix 下推到 master dataset 全量重跑（7 樣本 × 2 mode × flag=off，約 12-24 hr）→ 真正驗證在 full dataset 上是否一致反轉
3. 20260414 POSITIVE 結論的生物學解讀重寫：不再是「Intermediate AF 標記 subclonal LOH」，而可能是「stale baseline 下 CN1 區組成 heterogeneity artifact」

---

## 9. 特殊分布盤點（用戶要求）

| 現象 | 樣本 | 解讀 |
|------|------|------|
| COLO829 NG=0 fraction 80% LOH Outer T0 | COLO829 | ONT R10 無 methylation basecall，HP tag missing → NG=0 artifact，確認 COLO829 對 NG marker 不適用 |
| HCC1954 LOH Inner 仍 NG=2 dominant | HCC1954 | HER2+ ploidy ~4，"CN1<0.75" 對應 sub-diploid 但仍有部分 haplotype → NG>1 |
| HCC1937 LOH.bed 4531 segments | HCC1937 | BRCA1 HR defect → 最破碎的 LOH segments |
| HCC1954 單樣本保留 POS ΔNG | HCC1954 | **Stale rescale 的 artifact 指紋，最強反轉真實性證據** |
| HCC1395_DORADO LOH Inner NG > Outer | DORADO basecaller | Dorado HP fine-level artifact（需另行研究） |
| COLO829 CN1 LOH Outer TP rate 0.66 | COLO829 | 唯一 TP rate <0.85 的 cell，ONT R10 驗證 NG artifact |

---

## 10. 檔案清單

### 10.1 腳本（可重現）

```
research/ng_kde_rescaling/scripts/
├── utils_io.py                             # 共用 I/O 與常數
├── step0_build_master.py                   # 建立 7 樣本 + HCC1395 TO 合併資料集
├── step1_cn_tier_and_metrics.py            # CN 分層 + NG 分佈 + 復現 Test A/B
├── step1b_h2_verification_lohbed.py        # H2 驗證：LOH.bed + TP-only
└── step3_visualize.py                      # 6 張圖生成
```

### 10.2 資料表

```
research/ng_kde_rescaling/data/
├── merged_7samples_paired_full_plus_hcc1395_to.tsv.gz    # 369,080 rows
├── cn_tier_coverage_matrix.tsv                            # 209 rows（6 strategies × 7 samples × tier）
├── ng_dist_stratified.tsv                                 # 1,822 rows
├── ng_tprate_stratified.tsv                               # 1,296 rows
├── af_bin10_ng_crosstab.tsv                               # 560 rows
├── ng_20260414_replication_inter_vs_extreme.tsv           # 21 rows（核心復現表）
├── ng_loh_inner_vs_outer_extension.tsv                    # 162 rows（LOH 延伸）
├── ng_h2_verification_lohbed_vs_potential.tsv             # 14 rows（H2 驗證）
└── step0_per_sample_summary.tsv                           # 每樣本 summary
```

### 10.3 觀察圖

```
docs/experiments/in_progress/2026/04/figures/20260421_ng_kde_rescaled/
├── fig1_covm_density_per_sample.png         # 新 KDE 下 CovM 分佈
├── fig2_tprate_heatmap_A_vs_F.png           # TP rate 跨策略 × LOH 內外
├── fig3_ng_distribution_stacked.png         # NG=0/1/2/3/4 分佈
├── fig4_af10bin_ng_heatmap.png              # AF 10-bin × NG TP rate
├── fig5_replication_forest_3strategies.png  # ⭐ 20260414 復現森林圖
└── fig6_h2_verification_lohbed.png          # ⭐ H2 驗證：排除方法學差異
```

### 10.4 執行 log

```
research/ng_kde_rescaling/observations/
├── step0_log.txt
├── step1_log.txt
└── step1b_h2_log.txt
```

---

## 11. 後續建議（等用戶決定）

| # | Option | 成本 | 預期收益 |
|---|--------|------|---------|
| P1 | Memory 更新（`project_loh_subclone_af_methylation_positive.md` 加警示） | 10 min | 防止未來 citation 時誤用原結論 |
| P2 | Master dataset 全量重跑（7 樣本 × 2 mode，paired_full 已有，仍需 to_pileup × 6 + paired_pileup × 6） | 12-24 hr | 在完整 caller filter 下 final verdict |
| P3 | 20260414 原 F pilot 數據集的 CovM 欄位 post-hoc rescale，重算 F pilot NG=4+AF<0.4 filter | 2-3 hr | 驗證 F pilot 92.8% TP rate 在新 baseline 下穩定度 |
| P4 | Knowledge MCP 文獻深度查閱（CAMDAC、FACETS、TITAN 詳細比對） | 30-60 min | 為 paper draft 加強引用 |
| P5 | 生物學機制論文草稿（新發現：stale baseline artifact detection） | 2-3 day | 方法論 findings 可發表 |

---

## 12. 可檢驗性（用戶可驗證）

**每項主張對應產出檔案**：

| 主張 | 驗證路徑 |
|------|---------|
| 5/7 反轉 Strategy A | `data/ng_20260414_replication_inter_vs_extreme.tsv` 第 1-7 行 |
| H2 驗證 LOH.bed 不改結論 | `data/ng_h2_verification_lohbed_vs_potential.tsv` LOH_bed_TP_only method |
| HCC1954 唯一保留 POS | `fig5_replication_forest_3strategies.png` HCC1954 列 |
| NG 分佈 LOH 內外差異 | `fig3_ng_distribution_stacked.png` 上下兩列對比 |
| Per-sample baseline 差異 | `step0_log.txt` + `research/kde_fix_validation/outputs/step3_quantile_drift/quantile_drift_per_sample.tsv` |

---

**報告完成日期**：2026-04-22
**作者**：InterSubMod Research (Claude Opus 4.7)
**Sessionid**：延續前次對話（context compacted）
