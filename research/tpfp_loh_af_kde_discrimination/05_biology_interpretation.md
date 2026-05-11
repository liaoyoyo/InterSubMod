---
title: 生物學詮釋 — 為什麼這些 cell 富集 TP/FP
status: in_progress
last_updated: 2026-04-22
---

# 05. 生物學詮釋

本檔回答：**5D cube 中富集 TP 和富集 FP 的 cell，各代表什麼生物學情境？**
此詮釋**不做跨樣本泛化主張**（見 `06_limitations_future_work.md` §2 警示），僅解釋 HCC1395 TO 這個單一測試場的觀察。

---

## 5.1 高 TP rate cell 的生物學分類

將 Top-17 cells（見 `data/tpfp_5d_pareto_HCC1395_TO.tsv`）依生物情境分群：

### Group A — True LOH + Extreme AF（pure haplotype loss）

```
rank  LOH_Subtype    AF_class   cn_tier  NG  NR_band  n    TP%
 1    LOH_Strong     Extreme    T1       3   mid      22   100.0%
 3    LOH_Weak       Extreme    T2       3   mid      44    97.7%
 9    LOH_Strong     Extreme    T2       3   high     87    96.6%
15    LOH_Weak       Extreme    T1       3   mid      91    94.5%
```

**生物學解讀**：LOH 區的 variant，allele frequency 應該 ≈ 1.0（wild-type allele 已 lost）或 ≈ 0（variant allele 在 lost 側）。`AF_class = Extreme` + `LOH_Subtype` 非 None 的組合對應「純 haplotype」情境，FP 來源只剩：
- mapping artifact (低機率)
- reference variant pollution (極低機率)

**為什麼 NG=3**（而不是 NG=4）：NG=3 可能捕獲「somatic HP tag 在 tumor 上的三分群」— 即使 germline 只有兩 allele，HP tag 在 tumor region 可能因 somatic methylation 再分出一個 epigenetic sub-cluster。此即 `project_loh_subclone_af_methylation_positive.md` 記錄的現象。

**⚠ 警告**：NG≥3 訊號在 `--germline-hp-only` flag=on 下會消失（見 `project_readparser_germline_hp_only_phase1_negative.md`）。本觀察基於 flag=off 的值，其「real signal vs artifact」仍待 flag=on 重驗證。

### Group B — Non-LOH Diploid Het（canonical somatic het）

```
rank  LOH_Subtype  AF_class    cn_tier  NG  NR_band  n    TP%
 5    None         Near-half   T1       2   mid      171   97.7%
 8    None         Near-half   T2       3   high      31   96.8%
11    None         Near-half   T2       2   mid       21   95.2%
13    None         Near-half   T1       3   mid      114   94.7%
```

**生物學解讀**：這是 S3 Diploid Het 的核心。`AF ≈ 0.5` + `CN ≈ 2` + `None LOH` 對應 canonical somatic heterozygous variant — 一半 tumor reads 帶 variant allele，另一半帶 wild-type，基因組為 diploid。

**NG=2 vs NG=3 的細微差異**：
- NG=2（rank 5, 11）= variant/non-variant allele 各為一分群，NumReads 中等（60-120）
- NG=3（rank 8, 13）= 加入 somatic methylation 再分群，high NR 或 mid NR

這兩種情境都對應「純淨 somatic het」— 只是 NG tile 數因 methylation complexity 略異。

### Group C — LOH_Noise + Extreme + high_NR（意外發現）

```
rank  LOH_Subtype  AF_class  cn_tier  NG  NR_band  n    TP%
 2    LOH_Noise    Extreme   T3       2   high     56    98.2%
 6    LOH_Noise    Extreme   T2       2   high    146    97.3%
 7    LOH_Noise    Extreme   T2       3   high     35    97.1%
```

**這是最引人注目的新發現**：`LOH_Noise` 原被視為「LOH 訊號在 noise threshold 之下」的桶，應與 `None` 類似偏 FP 多。但這裡 TP rate 97%+。

**可能的生物學解讀**（需後續驗證）：
1. **深覆蓋揭示弱 LOH 訊號**：`nr_band = high` (≥120) 下，即使統計上被劃為 `LOH_Noise`，真實的 allele imbalance 仍存在 — 高覆蓋樣本的 noise threshold 保守。
2. **CN=3 (T3) amplified region with loss**：rank 2 的 T3 + Extreme 組合可能對應 `Amplified single-copy` — 單 copy amplify 後 AF 自然 extreme。此情境在 HCC1395 的 chr8 / chr1p 有先例（memory `project_hcc1395_chr8_hotspot.md`）。
3. **Sampling noise 被 high NR 克服**：低 NR 下 LOH 訊號會被歸為 Noise；same genomic event 在 high NR 樣本會升級為 Weak/Strong。所以 `LOH_Noise + high NR` 等同「被壓低分類的 real LOH」。

**實際驗證方式**（未做）：檢視這些 cells 落在 HCC1395 的哪些 chromosomes，是否集中在已知 LOH hotspots (chr3p, chr9p)。

### Group D — LOH_Weak + Intermediate + T1/T2（亞克隆訊號）

```
rank  LOH_Subtype  AF_class       cn_tier  NG  NR_band  n    TP%
 4    LOH_Weak     Intermediate   T1       2   mid      43   97.7%
10    LOH_Weak     Intermediate   T1       3   mid      27   96.3%
17    LOH_Weak     Intermediate   T2       2   high     24   91.7%
```

**生物學解讀**：`Intermediate AF` (0.1 < AF < 0.4 或 0.6 < AF < 0.9) + weak LOH 暗示亞克隆 (subclone) 訊號：
- 腫瘤非單一克隆，某個 subclone 佔 30-70% 並在該 region 有 LOH
- 這在 HCC1395 的高子克隆多樣性中已被 ISM 研究確認（見 `project_loh_subclone_af_methylation_positive.md`）

S2 Scheme（Subclonal_LOH_Inter）原本就設計捕捉此情境，但 Top-17 發現 S2 內**再細分 T1 subset**（CN=1-2 邊界）比整個 S2 更純：

| Cell | n | TP rate |
|------|--:|------:|
| S2 全體（LOH_Weak/Subclone + Intermediate）| 214 | 87.4% |
| Top-17 內 LOH_Weak+Intermediate+T1 subset | 70 (43+27) | 97.1% |

這說明 S2 本身還可以再切。

---

## 5.2 無判別力區塊（S4）的生物學解讀

S4 = None + Extreme AF，包含 30,432 variants，TP rate = 71.1%（= baseline）。這是整個 TO mode FP 的主要來源。

**生物學假設**：S4 的本質是 **germline leak** — TO mode 無 normal 參考，極端 AF (>0.9 or <0.1) 的 variant 可能是 germline homozygous 而被誤判為 somatic。

**證據**：
- fig_v2_9 (v2 §14) 顯示 S4 內 TP vs FP 的 AlleleDelta 分佈 **Cliff's δ = -0.49**（最強特徵差）
- AlleleDelta 在 FP 偏更極端 → 對應 germline-only methylation pattern（normal 與 tumor 一致）
- 單特徵 threshold (E5) 最多 fold 1.17× — 不足以過濾

![v2_9 S4 內 TP vs FP feature violin](figures/existing/fig_v2_9_feature_violin_in_S4.png)

**這是為什麼 S4 無法 rescue**：germline leak 在 TO mode 天生無法區分，除非引入 normal 參考（即回到 paired mode）。這回到本專案的更大定位：**ISM 的價值在 characterization，不在 variant filter**。

---

## 5.3 Top-17 在 chromosome 上的分佈（obs04）

- Top-17 cells 分散於 22 條 autosome + chrX
- 無明顯 chromosome bias（與 HCC1395 全域 LOH 分佈一致）
- 這**排除**「Top-17 是單一 chromosome artifact」的假設

若 Top-17 集中在 chr8 或 chr3p（HCC1395 known LOH hotspots），會擔心 biology 被 LOH hotspot 驅動；實際上散佈於所有染色體。

![obs04 Top-17 / Top-28 / S3 在 chromosomes 上的空間分佈](figures/new/obs04_chr_spatial_scheme.png)

![v2_6 生物模組詮釋 heatmap](figures/existing/fig_v2_6_biology_module_interpretation.png)

---

## 5.4 為什麼 S3 vs Top-17 的差異重要

S3 (380 regions) 和 Top-17 (1099 regions) 的 337 個共同 region 是「canonical somatic het」群；多出的 762 個 Top-17-only regions 屬於：

- Group A（LOH + Extreme + NG≥3）= ~250 regions
- Group C（LOH_Noise + Extreme + high NR）= ~230 regions
- Group D（LOH_Weak + Intermediate + T1）= ~100 regions
- 其他 zero-FP cells（n≥20 但 n_fp=0）= ~180 regions

**生物學意涵**：Top-17 代表「**S3 canonical het + LOH-touched extreme/intermediate 區**」的聯集。S3 只切到最保守的 canonical 部分；Top-17 把 LOH 相關的高純度區一併納入。

---

## 5.5 COLO829 座標不同的生物學解讀

COLO829 paired_full top cells 主要在 NG=1 + Intermediate AF。為什麼差異這麼大？

| 因素 | HCC1395 TO | COLO829 paired_full |
|------|-----------|--------------------|
| Baseline TP% | 71.1% | 93.9% |
| FP 主要來源 | Germline leak (S4) | LOH boundary + CN-neutral LOH |
| 高純度區主軸 | NG=2-3（methylation complexity）| NG=1（single haplotype simple）|
| AF 主軸 | Near-half + Extreme | Intermediate + Extreme |
| Tumor 特性 | Triple-negative breast ca, 高子克隆 | Melanoma, 較簡單克隆結構 |

**推測**：COLO829 的 FP 熱區不在 S4（因為 baseline 已 94%），而在 subtle LOH 邊界；其 top cells 捕獲的是 **LOH boundary-adjacent Intermediate AF**，需要 NG=1 （單一 haplotype HP tag）作為乾淨訊號。

HCC1395 TO 則是 germline leak 大量，top cells 需要透過 **NG=2-3 多 haplotype 訊號** 證明是真 somatic event。

**這兩個樣本的 biology 結構不同**，top cells 座標自然不同 — 但「**5D cube 存在 Pareto-optimal 點**」這個 meta-claim 在兩個樣本都成立。

![obs06 per sample NG 分佈 TP vs FP（顯示 HCC1395 TO vs COLO829 NG 軸差異）](figures/new/obs06_ng_by_sample_stacked.png)

![v2_11b COLO829 paired_full 5D panorama](figures/existing/fig_v2_11b_panorama_COLO829.png)

---

## 5.6 關鍵生物學假說（需進一步驗證）

| 假說 | 現況 | 驗證方式 |
|------|-----|---------|
| Group C (LOH_Noise + high NR) 是真實 LOH | 純觀察 | 在 `tp_pileup` 純 pileup 資料下驗證 LOH 訊號連續性 |
| Group A 的 NG=3 是 methylation-confirmed | 依賴 flag=off 值 | flag=on 重跑看 NG=3 是否消失 |
| S4 無法 rescue 是 germline leak | AlleleDelta 證據（fig_v2_9）| 在 paired mode 比對 normal read 直接確認 |
| Top-17 vs S3 的多餘 762 regions 是 LOH-adjacent | obs04 空間分佈分散 | 疊加 LOH.bed 座標驗證 |

---

## 5.7 與論文定位的關聯

ISM 在 2026-04 確認定位為「**read-level epigenetic characterization**，非 variant filter」（見 CLAUDE.md「開發重點」）。本研究的生物學詮釋與此定位**一致**：

- Top-17 並非最優 filter（recall 只有 3.7%，無法取代 caller）
- Top-17 的價值在於**標記高純度 somatic event 的 epigenetic 情境**（Group A-D 四類 biology）
- 對於 variant filter 任務，Top-17 僅能作為**高純度白名單 boost**，不能作為 precision 主力

未來工作：將 Top-17 的四個 Group 轉化為**臨床相關的 epigenetic signatures**（如 LOH-associated methylation pattern 白名單），而非單純的 TP/FP filter。
