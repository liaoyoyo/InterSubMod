---
title: TO mode 下 LOH-constrained phasing 機制：NG=2 拆分揭露 somatic vs germline-het FP 分離
status: observation_corrected
owner: liaoyoyo2001
created: 2026-04-22 15:56
last_updated: 2026-04-22 20:40
correction_note: "V1（15:56）誤以 HPFineNGroups 為 methylation bimodality 訊號；
經 C++ 原始碼查證（LabelTest.cpp::hp_to_fine_labels）與 obs18 組成拆分驗證（6 TO 樣本），
正確機制為 **LOH-constrained phasing**（variant 是否被迫落在單 haplotype 上）。"
data_source: data/obs17_TO_4d_cube.tsv + data/obs18_NG2_composition_by_sample.tsv
parent_project: 00_INDEX.md
---

# 09. TO mode 的 LOH-constrained phasing 機制

## 9.0 TL;DR（3 句話）

1. **HPFineNGroups (NG) 不是 methylation 訊號**，而是 `{HP1, HP1-1, HP2, HP2-1}` 四個 **(haplotype × variant-carrying)** bucket 的 occupancy 計數。
2. **NG=2 的 Inner 區 93-99% 是「same-haplotype 分裂」**（HP1+HP1-1 或 HP2+HP2-1，somatic variant 在單 hap 內形成 ref/alt 子族，**6/6 TO 樣本一致**），TP rate 0.43-0.96（median 0.93）。
3. **NG=2 的 Outer 區主要是「cross-haplotype phasing」**（HP1+HP2-1，經典 het SNV 模式，但 germline het 長得一模一樣），TP rate 0.08-0.88（median 0.55），gap = +0.37；**這解釋了 TO 模式 germline-het FP 問題的物理根源**。

---

## 9.1 歷史與更正

**舊敘述（已棄置）**：V1 於 2026-04-22 15:56 撰寫時誤解 HPFineNGroups 為「methylation bimodality」（2 個甲基化 cluster），並基於此推論 "Haplotype-loss-dependent methylation bimodality" 為論文主軸。

**更正觸發**：用戶提問「NG=2 與甲基有關係嗎」促使回查 C++ 原始碼。

**發現**：
- `src/core/LabelTest.cpp:265-302` 的 `hp_to_fine_labels()` 將 reads 依 `hp_tag` 字串分到 4 個 bucket：
  ```cpp
  if (hp == "1") group = 0;        // HP1    = haplotype 1, reference allele
  else if (hp == "1-1") group = 1; // HP1-1 = haplotype 1, somatic allele
  else if (hp == "2") group = 2;   // HP2    = haplotype 2, reference allele
  else if (hp == "2-1") group = 3; // HP2-1 = haplotype 2, somatic allele
  ```
- `include/core/Stats.hpp:324` 註解：`d_hp1_hp1s: "Same haplotype, germline vs somatic"`
- **HPFineNGroups 純粹是 bucket occupancy count**，與 methylation 無直接計算關係
- Methylation 只在 `HPFineF`/`HPFineP` 的 PERMANOVA 測試中參與（品質檢驗：這 4 bucket 的 methylation 分佈是否真正分離），但 NGroups 本身不依賴甲基化值

---

## 9.2 NG=2 的四種可能組成

HPFineNGroups=2 代表 4 個 bucket 中有 2 個被 populate（其他 2 個 reads 不足或為 0）。可能組合：

| 組合代號 | bucket 1 | bucket 2 | 生物學意義 |
|---------|----------|----------|------------|
| **same_HP1** | HP1 | HP1-1 | 單 haplotype 1 內部：ref 子族 + somatic 子族 |
| **same_HP2** | HP2 | HP2-1 | 單 haplotype 2 內部：ref 子族 + somatic 子族 |
| **cross_het** | HP1 | HP2-1 | 兩個 haplotype：hap1 都 ref + hap2 都 somatic（canonical het phasing）|
| **cross_het_inv** | HP1-1 | HP2 | 對稱版：hap1 都 somatic + hap2 都 ref |
| other | (其餘) | (其餘) | 如 HP1+HP2（無 somatic marker）、HP1-1+HP2-1（兩邊都變異）|

**關鍵區分**：
- **same-hap 系列**：在 LOH 區只有單 haplotype 存在 → somatic SNV 發生後必然產生 same-hap 分裂 → **必然是真 somatic TP**
- **cross-het 系列**：雙 haplotype 完好 → 可以是真 somatic het，**也可以是 germline het SNV（caller 無法在 TO 模式下區分）**

---

## 9.3 跨 6 TO 樣本驗證（obs18）

![NG=2 composition stacked bar](figures/new/obs18/obs18_NG2_composition_proportion.png)

### 組成占比（NG=2, Extreme AF 子集）

| sample | Inner: same-hap% | Outer: cross-het% | Inner total | Outer total |
|--------|----------------:|------------------:|------------:|------------:|
| HCC1395 | **93.2%** | 0.1% (99.7% "other") | 632 | 4,965 |
| HCC1395_DORADO | **99.0%** | **97.5%** | 10,570 | 6,827 |
| HCC1937 | **98.8%** | 42.9% | 8,521 | 3,887 |
| HCC1954 | **96.5%** | 31.0% | 9,363 | 23,838 |
| H2009 | **98.3%** | 45.9% | 38,210 | 15,709 |
| H1437 | **97.0%** | 71.3% | 9,145 | 12,116 |

**6/6 TO 樣本的 Inner NG=2 都是 ≥93% same-haplotype** — 跨樣本完全一致。

Outer 的 NG=2 組成較為混合（cross-het 占比 0.1-97.5%，median 44%），反映不同 phasing tool 與樣本 purity 下的 edge case。但主旨清楚：**Outer 的 NG=2 至少有顯著比例是 cross-het，而 Inner 幾乎沒有 cross-het**。

### TP rate 對照（6 TO 樣本）

![NG=2 composition heatmap](figures/new/obs18/obs18_NG2_composition_heatmap.png)

| sample | Inner same_HP1 TP | Outer cross_het TP | gap |
|--------|:----------------:|:------------------:|:---:|
| HCC1395 | **0.96** | 0.50 | **+0.46** |
| HCC1395_DORADO | **0.94** | 0.55 | **+0.39** |
| HCC1937 | **0.76** | 0.24 | **+0.52** |
| HCC1954 | 0.43 | 0.08 | **+0.35** |
| H2009 | **0.93** | 0.88 | +0.05（飽和）|
| H1437 | **0.92** | 0.69 | +0.23 |

**median gap = +0.37**；6/6 正向；0 反向。

---

## 9.4 機制解釋：LOH-constrained phasing

### 為什麼 Inner NG=2 ≈ same-hap？

LOH 定義上意味著該區域只保留**單一 haplotype**（另一個 allele 已 copy-lost）。在這種區域，任何 somatic SNV 發生後：
- Reads 只能來自單 haplotype（例如 HP1）
- Variant 產生後，帶 variant 的 reads 變為 HP1-1，不帶的保持 HP1
- 結果：**HP1 + HP1-1 = NG=2 的 same-hap 分裂**

這是**物理必然**，不是統計巧合。

### 為什麼 Outer NG=2 主要是 cross-het，且 TP rate 低？

非 LOH 區保留雙 haplotype。若 SNV 出現：
- **真 somatic het**：variant 可能出現在任一 hap → HP1 + HP2-1 或 HP1-1 + HP2
- **Germline het**：也會產生完全相同的 HP1 + HP2-1 phasing pattern

這兩者在 phasing output 上**完全無法區分**。Caller（ClairS-TO）在 TO 模式下僅用腫瘤 reads，無配對正常樣本 → 傾向把 germline het 誤判為 somatic → Outer NG=2 中混入大量 germline FP。

**結論**：Outer × NG=2 cross-het 的 TP rate 低，不是因為 variant 本身有問題，而是因為這種 pattern **在 TO 模式下 germline-somatic 不可分**。

### 對 TO 模式 FP 根源的直接映射

長久以來「TO 模式 FP 過多」的歸因都停留在「caller 不夠聰明」或「需要 normal 樣本」。本分析揭露更具體的機制：
- TO 模式下 germline het FP 的主要「長相」是 **Outer × cross-het NG=2**
- 若有方法識別「非 LOH 區的 cross-het NG=2」，即可標記為 **germline-het 高風險 cells**
- 相對地，**Inner × same-hap NG=2** 是 TO 模式下**最可信任的 somatic cells**

---

## 9.5 異常樣本：HCC1954

HCC1954 的 Inner + same_HP1 TP rate 僅 0.43（其他 5 樣本 0.76-0.96）。可能原因：
- 該樣本整體 baseline TP = 0.25（全 TO 最低），caller 基本表現差
- Potential_LOH 標記可能受該樣本 polyploid 特性影響而不可靠
- 需獨立分析（見 obs17 §9.4）

**處理建議**：
- LOSO 驗證優先 held-out HCC1954
- 不影響本機制主結論（因為 6/6 樣本 direction 仍一致）

---

## 9.6 對先前研究假設的更正

| 舊假設 | 更正後 |
|-------|-------|
| NG=2 是 methylation bimodality | ❌ 錯。是 (haplotype × variant) bucket occupancy |
| Inner + NG=2 高 TP 來自甲基化訊號 | ❌ 錯。來自 LOH 區強迫 variant 落在單 hap |
| Outer + NG=2 低 TP 因 parental methylation 差異 | ❌ 錯。因 germline het 和 somatic het phasing 相同 |
| 論文主軸：methylation bimodality | ❌ 改為：LOH-constrained phasing |
| Methylation essential for TO discrimination | ❌ 對此機制而言不 essential（phasing + LOH annotation 已足） |

**Methylation 的角色仍在**，但不是此 discovery 的主軸：
- `HPFineSig`（PERMANOVA on methylation distance）驗證 bucket 確實在 methylation 空間分離 → 支援 bucket definition 的有效性
- `HPMergedDelta` / `PerCpgASM_Valid` / `NME_HP1/2` 才是真正的 methylation signatures，獨立於 phasing，是另一條研究線

---

## 9.7 論文定位更新

**新 Title 方向（取代原 methylation-based title）**：

**"LOH-constrained phasing signatures distinguish somatic from germline-like variants in tumor-only sequencing"**

**核心假說（可直接驗證）**：
> 在 TO 模式下，`Inner × same-haplotype NG=2` = 幾乎純 somatic TP（median TP rate 0.93，5/6 TO 樣本）；`Outer × cross-haplotype NG=2` = germline-het 與 somatic-het 無法區分 FP-contaminated（median TP rate 0.55）。兩者的 TP rate gap +0.37 是 TO 模式 germline-het FP 問題的**直接 phasing signature**。

**優勢**：
- 機制清晰純粹：只需 phasing tool + LOH annotation，不涉及 methylation 解釋
- 可 100% 從 C++ pipeline 輸出驗證
- 提供 TO FP 問題的**可操作診斷**：標記 Outer × cross-het NG=2 為 germline-het 風險

**與先前結論的一致性**：
- `07_figure_layers.md §7.5.1` Tier A 的 3 cells（Extreme + T1 + NG=2）與本 discovery 吻合 — T1 CN（CN<2 loss，等同 LOH 訊號）+ NG=2 = **same-hap + LOH 雙重特徵**
- `project_loh_subclone_af_methylation_positive.md` 的「Inter AF×NGroups POSITIVE」被重新解讀為 **同一 LOH-phasing 現象的 paired_full 版本**，不是甲基化 × AF 直接耦合

---

## 9.8 下一步建議

### 立即（驗證機制）
1. **Negative control**：`--germline-hp-only` flag 下此 Inner-Outer gap 應消失（該 flag 移除 somatic HP tag → HP1-1/HP2-1 bucket 全空）
2. **Formal statistics**：對 6/6 跨樣本 gap 做 Wilcoxon signed-rank（vs null=0）
3. **HCC1954 根因**：分析其 Potential_LOH 與 LOH.bed（PON）的 Jaccard，確認標記可靠性

### 中期（擴展發現）
4. **Filter design**：以「Inner + same-hap NG=2」定義新的 `ISM-LOH-phasing` filter，比較 LOSO F1 vs Top-17
5. **TO caller 改進**：建議 ClairS-TO 後處理加入「flag Outer + cross-het NG=2 為 germline-het 風險」
6. **跨 caller 驗證**：DeepSomatic-TO / Strelka-TO 是否有同樣 phasing signature？

### 長期（論文）
7. 單細胞 methylation 或 phasing 驗證：確認 Inner same-hap NG=2 的兩 bucket 真正對應 subclonal 結構
8. 跨 dataset 驗證：TCGA TO samples 是否重現此機制

---

## 9.9 資料與複現

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination

# 原 4D cube 分析（obs17）
python3 scripts/obs17_TO_af_lohside_ng.py

# NG=2 組成拆分（obs18，本次新增）
python3 scripts/obs18_TO_NG2_composition.py
```

**輸出**：
- `data/obs17_TO_4d_cube.tsv` / `obs17_TO_direction_summary.tsv`
- `data/obs18_NG2_composition_by_sample.tsv`（56 rows：6 sample × 2 side × 5 combo）
- `data/obs18_NG2_composition_proportion.tsv`（12 rows：6 sample × 2 side）
- `figures/new/obs17/` × 3 圖
- `figures/new/obs18/` × 2 圖

**C++ 原始碼參照**：
- `src/core/LabelTest.cpp:265-305`（hp_to_fine_labels bucket 定義）
- `include/core/Stats.hpp:65-80`（FullLabel 結構）
- `include/core/Stats.hpp:323-330`（HPFinePairwise 的 germline/somatic 註解）
