<!--
建立時間: 2026-04-02 16:00
目標: 以 HCC1395 purity 梯度 subsample 驗證 self-phasing circular dependency 在不同腫瘤純度下的表現
處理範圍: HCC1395 subsample 5 purities (20%-100%) × 22 chromosomes × SEQC2 truth set
關聯檔案:
  - docs/reports/validated/2026/04/20260402_longphase_to_vs_s_causal_chain_report_01.md
  - scripts/analysis/build_purity_self_phasing_vcf_analysis.py
  - scripts/analysis/build_purity_loh_proxy_analysis.py
  - scripts/analysis/plot_purity_self_phasing_figures.py
-->

# Self-Phasing Circular Dependency 的 Purity 依賴性驗證

## 動機

先前因果鏈報告（`20260402_longphase_to_vs_s_causal_chain_report_01.md`）在 **100% purity 癌症細胞株**上確認了 self-phasing circular dependency 是 TO-only LOH 的主因。但報告中 Limitation 6.3.1 指出：

> 臨床樣本的 tumor purity 通常為 30-70%，LOH enrichment 翻轉的幅度可能顯著減弱，因為低 purity 時 somatic allele 的 HP 偏移效應被 normal reads 稀釋。

本報告利用 HCC1395 已有的 **5 個 purity 梯度 subsample**（20% ~ 100%）的 ClairS-TO phased VCF 輸出，直接驗證此假設。

---

## 資料來源

| Purity | Subsample | 估計 VAF | ClairS-TO PASS | SEQC2 TP |
|:------:|:---------:|:--------:|:--------------:|:--------:|
| 100% | t50_n00 | 0.50 | 1,831,636 | 27,299 |
| 80% | t40_n10 | 0.40 | 2,448,704 | 28,766 |
| 60% | t30_n20 | 0.30 | 2,675,135 | 26,648 |
| 40% | t20_n30 | 0.20 | 2,702,531 | 21,233 |
| 20% | t10_n40 | 0.10 | 2,668,849 | 9,734 |

- **Truth set**: SEQC2 high-confidence sSNV v1.2.1（39,447 positions）
- **Phased VCF**: ClairS-TO v0.3.0 內建 LongPhase-TO phasing 輸出（per-chromosome）
- **分析範圍**: chr1-chr22（autosome only）

---

## 結果一：Self-Phasing Scaffold 侵入率隨 Purity 線性下降

Self-phasing 的前提是 somatic variant 進入 phasing scaffold（有 PS 值）。隨著 purity 下降，somatic VAF 降低，其 ALT reads 對 phasing graph 的貢獻減弱，導致進入 scaffold 的比例下降。

![Figure 1: TP Scaffold Rate vs Purity](purity_figures/fig1_tp_scaffold_rate_vs_purity.png)

| Purity | TP 在 Scaffold 中 | 含義 |
|:------:|:-----------------:|:-----|
| 100% | **60.2%** | 超過半數 TP somatic 作為 phasing anchor |
| 80% | 56.2% | 輕微下降 |
| 60% | 45.0% | 顯著下降至一半以下 |
| 40% | 32.5% | 僅三分之一 |
| 20% | **28.9%** | 不足三成 |

**趨勢**：scaffold 侵入率從 60% → 29%，降幅超過一半。低 purity 下 somatic variant 對 phasing 的汙染顯著減輕。

---

## 結果二：Scaffold LOH 效應在 ≤60% Purity 消失

核心問題：進入 scaffold 的 TP 是否比未進入 scaffold 的更容易出現 LOH-like（AF > 0.9）？

![Figure 2: Scaffold LOH Effect vs Purity](purity_figures/fig2_scaffold_loh_effect_vs_purity.png)

| Purity | Scaffold TP LOH% | Non-scaffold TP LOH% | 放大倍數 | 嚴重度 |
|:------:|:-----------------:|:--------------------:|:--------:|:------:|
| 100% | **16.61%** | 0.64% | **26.1×** | 嚴重 |
| 80% | 0.51% | 0.06% | **8.0×** | 可測量 |
| 60% | 0.008% | 0.007% | 1.2× | **可忽略** |
| 40% | 0% | 0% | — | 無 |
| 20% | 0.04% | 0% | — | 無 |

**關鍵發現**：

- **100% purity**：scaffold 內 TP 的 LOH-like 率是外部的 **26 倍** — self-phasing 造成嚴重 LOH artifact
- **80% purity**：效應仍在（8 倍），但絕對值已降至 0.5%
- **≤60% purity**：放大倍數降至 1.2 倍 — **self-phasing LOH 效應消失**

---

## 結果三：TP AF 分佈確認 — 低 Purity 無法達到 LOH-like 閾值

Self-phasing LOH 需要 somatic variant 有足夠高的 AF。隨 purity 下降，TP AF 被稀釋，能達到 AF > 0.8 的 TP 急劇減少。

![Figure 3: TP AF Distribution by Purity](purity_figures/fig3_tp_af_distribution_by_purity.png)

| Purity | TP AF median | 理論 VAF (P/2) | TP AF>0.8 數量 | TP AF>0.8 比例 |
|:------:|:-----------:|:---------:|:---------:|:---------:|
| 100% | 0.412 | 0.50 | 3,409 | **12.5%** |
| 80% | 0.354 | 0.40 | 1,165 | 4.0% |
| 60% | 0.279 | 0.30 | 27 | 0.1% |
| 40% | 0.214 | 0.20 | 5 | **0.0%** |
| 20% | 0.154 | 0.10 | 2 | **0.0%** |

**解讀**：在 ≤60% purity 下，幾乎沒有 TP somatic variant 能達到 AF > 0.8，因此即使它們進入了 scaffold，也不會產生 LOH-like 信號。Self-phasing 的「燃料」（高 AF somatic reads）被 normal reads 稀釋殆盡。

---

## 結果四：VCF AF vs ISM HP_Ratio — 兩種指標方向相反

VCF-level AF 分析顯示 LOH-like（AF > 0.9）在**所有 purity 都是 FP-enriched**。這與因果鏈報告中 ISM HP_Ratio 分析的 **TP-enriched** 方向相反。

![Figure 4: LOH Enrichment Direction vs Purity](purity_figures/fig4_loh_enrichment_direction_vs_purity.png)

| Purity | TP LOH@0.9 | FP LOH@0.9 | FP/TP 倍數 | 方向 |
|:------:|:---------:|:---------:|:--------:|:----:|
| 100% | 10.26% | 41.03% | **4×** | FP-enriched |
| 80% | 0.31% | 26.48% | **85×** | FP-enriched |
| 60% | 0.008% | 17.92% | **2388×** | FP-enriched |
| 40% | 0% | 16.42% | — | FP-enriched |
| 20% | 0.01% | 15.65% | **1524×** | FP-enriched |

**這個「矛盾」恰好證明 self-phasing 的機制**：

- **VCF AF** 測量 allele frequency（REF vs ALT reads 比例）→ FP 很多是 germline homozygous（AF ≈ 1.0）→ 天然 FP-enriched
- **ISM HP_Ratio** 測量 haplotype balance（HP1 vs HP2 reads 比例）→ self-phasing 讓 TP ALT reads 偏向一個 HP → **TP-enriched**

兩個指標方向相反 = **self-phasing 作用於 HP tag assignment 而非 allele frequency 本身**。

---

## 結果五：Phase Block 結構隨 Purity 改變

![Figure 5: Phase Block Structure vs Purity](purity_figures/fig5_phase_block_structure_vs_purity.png)

| Purity | Phase Blocks 數 | Median Size | Block 被 TP 汙染率 |
|:------:|:---------------:|:-----------:|:------------------:|
| 100% | **4,072** | **80** | 60.6% |
| 80% | 1,590 | 498 | 72.5% |
| 60% | 1,257 | 667 | 69.1% |
| 40% | 1,224 | 727 | 65.4% |
| 20% | 1,286 | 660 | 52.4% |

**意外發現**：100% purity 有**最多 blocks 但最小 median size**。原因：純腫瘤中大量 somatic variants 進入 scaffold → phasing graph 更複雜 → 產生更多但更短的 phase blocks（phasing fragmentation）。低 purity 時 germline 主導 scaffold → 更少、更穩定的 blocks。

---

## 綜合總結圖

![Figure 6: Combined Summary](purity_figures/fig6_combined_purity_summary.png)

---

## 結論

### Self-phasing 效應的 Purity 依賴性 — 三段分級

| Purity 範圍 | Self-phasing 嚴重度 | Scaffold LOH 放大 | 建議 |
|:-----------:|:------------------:|:-----------------:|:----:|
| **≥90%**（cell line / 高純度腫瘤） | **嚴重** | 26× | 必須校正 |
| **70-90%**（高純度臨床） | **可測量** | ~8× | 建議校正 |
| **≤60%**（典型臨床 purity） | **可忽略** | ≤1.2× | 無需特別處理 |

### 核心機制

```
Self-phasing 效應強度 ∝ somatic VAF ∝ tumor purity

  Purity ↓ → VAF = P/2 ↓ → ALT reads 不足以主導 phasing anchor
    → Scaffold 侵入率 ↓（60% → 29%）
    → LOH-like artifact ↓（16.6% → 0%）
    → Self-phasing 問題消失
```

### 對因果鏈報告 Limitation 的回應

> **原假設**：「臨床樣本的 tumor purity 通常為 30-70%，LOH enrichment 翻轉的幅度可能顯著減弱」

**驗證結論**：**假設成立且效應比預期更強** — 在 ≤60% purity 下 self-phasing LOH 效應不只是「減弱」，而是**基本消失**。我們的 100% cell line 結論代表 **worst case scenario**。

但 TO 模式的其他固有問題（31% germline FP 汙染、缺乏 normal baseline、HP tag 品質差異）在各 purity 下仍然存在。

---

## 數據與腳本

### 分析腳本
- `scripts/analysis/build_purity_self_phasing_vcf_analysis.py` — Scaffold 侵入率 + AF 分層
- `scripts/analysis/build_purity_loh_proxy_analysis.py` — LOH proxy + enrichment 分析
- `scripts/analysis/plot_purity_self_phasing_figures.py` — 6 張圖表

### 輸出數據（workspace: `20260402_phasing_causal_chain/`）
| 檔案 | 內容 |
|------|------|
| `purity_self_phasing_vcf_summary.tsv` | 5 purities 主要指標彙總 |
| `purity_self_phasing_af_stratified.tsv` | AF 分層 × purity 交叉表 |
| `purity_loh_proxy_by_af_threshold.tsv` | LOH proxy TP/FP enrichment |

### 圖表（`purity_figures/`）
| 圖表 | 內容 |
|------|------|
| `fig1_tp_scaffold_rate_vs_purity.png` | TP scaffold 侵入率 vs purity |
| `fig2_scaffold_loh_effect_vs_purity.png` | Scaffold LOH 放大效應 + ratio |
| `fig3_tp_af_distribution_by_purity.png` | TP AF 分佈堆疊圖 |
| `fig4_loh_enrichment_direction_vs_purity.png` | VCF AF LOH enrichment 方向 |
| `fig5_phase_block_structure_vs_purity.png` | Phase block 數量 + median size |
| `fig6_combined_purity_summary.png` | 四合一綜合總結圖 |
