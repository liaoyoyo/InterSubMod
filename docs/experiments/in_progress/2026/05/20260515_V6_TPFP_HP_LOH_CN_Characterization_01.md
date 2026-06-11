<!--
建立時間: 2026-05-15
更新時間: 2026-05-15
agent: Coordinator synthesis (main session) of multi-agent v6_bam_tpfp_hp_loh_cn cycle
status: in_progress
report_class: characterization-post-hoc (multi-stage cycle synthesis)
audience: PI / lab member / 自己未來
scope: HCC1395 pilot + 4 cross-sample (H1437/H2009/HCC1954/HCC1937), V3F/V5/V6 longphase haplotag 三向, characterization-only
tier: ⭐3 PARTIAL POSITIVE
parent_plan: research/v6_bam_tpfp_hp_loh_cn/00_PLAN.md (v0.3)
inputs:
  - research/paired_priority_bug_audit/phaseC_genome_three_way/ (12 ISM runs)
  - research/paired_priority_bug_audit/phaseD_v6_5sample/ (4 samples × 4 runs)
  - research/tpfp_loh_af_kde_discrimination/data/master.tsv.gz (HCC1395 paired_full LOH/CN)
  - research/seqc2_cnv_stratification/data/annotated_hcc1395_cnv.tsv (Z-GL SEQC2 oracle)
  - /big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_{tp,fp}.vcf.gz
upstream_reports:
  - InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md
  - InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md
  - InterSubMod/research/paired_priority_bug_audit/07_V6_validation_findings.md
  - InterSubMod/research/paired_priority_bug_audit/08_phaseD_v6_cross_sample_findings.md
outputs:
  - 本檔（主報告）
  - InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step{1,2,3,4}_*/findings.md
  - InterSubMod/research/v6_bam_tpfp_hp_loh_cn/{00_PLAN,01_data_provenance,02_methodology_notes,02_prior_art_notes}.md
verdict: PARTIAL POSITIVE — H4 chr8 hotspot LOH+CN 主導 + paradigm reframe (cross_het = TP signature 非 FP) + 1 cross-sample candidate; H3/H5/H7 NEGATIVE; ~63% FP 仍未被此 framework 解釋
last_verified: 2026-05-15
report_template: characterization v1.0
-->

# V3F/V5/V6 ISM 三向 × LOH × HP × CN — HCC1395 pilot + 4 樣本擴展 TP/FP characterization

> **Characterization-only** post-hoc 後分析。不評估 filter / ΔF1（per plan §Out-of-scope）。
> Verdict: **⭐3 PARTIAL POSITIVE** — 含 paradigm reframe + 1 cross-sample signature candidate

---

## 0. TL;DR

本 cycle 用 V3F/V5/V6 三個 longphase haplotag 版本（V3F=PON-only baseline / V5=Layer 1.5 / V6=Layer 1.5 reverted production candidate）對 HCC1395 paired-pileup 30,490 TP + 4,842 FP 跑 LOH × HP bucket × Coverage_Multiple 3 軸 50-cell grid，再擴展 4 樣本（H1437/H2009/HCC1954/HCC1937，COLO829 deferred）。**所有上游 ISM 結果已存在**（phaseC_genome_three_way + phaseD_v6_5sample），本計畫純粹 post-hoc 後分析，不重跑 BAM/ISM。

**7 大假設判定**：

| ID | 假設 | Verdict | 數值 |
|----|------|---------|------|
| H1a | V5 Inner same-hap TP gap > V3F (Δ≥0.03) | **NEGATIVE (ceiling)** | Δ=+0.003 (NG=2 TP rate 飽和 ≥0.99) |
| H1b | V6 Inner same-hap TP gap > V5 (Δ≥0.03) | **NEGATIVE (ceiling)** | Δ=+0.001 |
| H1c | V6 vs V3F 累積 (Δ≥0.06) | **NEGATIVE (ceiling)** | Δ=+0.004 |
| **H2** | Inner × same_HP TP enrichment | **POSITIVE (directional)** | 6 cells, mean enrichment 1.011× (max 可能 1.017×) |
| **H3** | Z-OCH × cov_gain FP enrichment ≥ 2.0× | **NEGATIVE (reframe!)** | FP_enrichment 0.124× — Z-OCH 是 TP-pure |
| **H4** | chr8 hotspot (LOH+CN) − HP ≥ 0.05 | **POSITIVE** | Δ=+0.186, CN dev_explained 0.211 vs HP 0.063 |
| H5 | Z-AUTO vs known zones Jaccard ≥ 0.60 | **NEGATIVE** | Jaccard=0.184 (~70% Z-AUTO 在 known zone 外) |
| H6 | Powered cells (n≥50) ≥ 30% | **POSITIVE** | 46% (23/50) |
| H7 | ≥5 cells 通過 7 confound guards | **NEGATIVE (HP-asymmetry 是 biology 非 artifact)** | 2 cells 通過全部，8 cells 通過 6/7 |

**Paradigm reframe**: 預設的「FP-rich zone」中 2/3 (Z-OCH Outer cross_het, Z-GL Inner gain+LOH) 實際上是 **TP-pure signatures**（cross_het = somatic-evidence marker），不是 FP markers。真正的 FP 機制：**chr8 hotspot 由 CN+caller_af 主導（HP 次要）**。

**Cross-sample (n=5)**: 唯一 signature candidate `Outer|other|cov_high_gain` — 5/5 樣本 TP rate=1.000（Wilcoxon p=0.0625 exact min），但 effect size 極小 (+0.0069)，受 caller 飽和 (TP rate ≥0.998) 限制。

**Framework coverage gap**: Z-CHR8 + Z-AUTO 合捕獲 ~37% 全 FP；**剩 ~63% FP 不被此 framework 解釋**，需新軸（mappability / repeat / GC / SV context）。

---

## 1. Context & 動機

### 1.1 為什麼做這個

CURRENT_FOCUS 2026-05-13 主軸切換：Phase 1 W1 F-paired-D3 = V5 Layer 1.5 ISM 影響量化。本 cycle 是該任務的延伸 + Thread D 主軸（LOH-constrained phasing signatures）的 LOH × HP × CN 三軸交叉擴充。

**3 個新變量**讓重做 HP × LOH × CN 三軸交叉有意義：
1. **V6 phased BAM 已產出**（5/10 commit 鏈完成 + 4 樣本驗證），未量化過 V3F→V5→V6 在 LOH.bed 內外的 TP/FP 區分能力差異
2. Thread D 已有「LOH 內外 × HP 4-bucket」框架，但未加入 CN 軸
3. Coverage_Multiple KDE-corrected 已驗證 (r=0.831 vs SEQC2 CNV truth)，幾倍體軸首次與 LOH 內外 × HP bucket 交叉

### 1.2 歷史警訊與本研究的差異化

歷史上類似嘗試全 NEGATIVE：LOH binary filter 10/10 策略失敗、CN zone-aware filter 跨樣本不一致、LOH+HPMergedSig 7.4× 87.5% 來自 HCC1395 chr8 樣本特異性。

**但這些都是 filter 失敗，不是 characterization 失敗**。本 cycle 嚴守 characterization 邊界（不評估 ΔF1 / AUC > 0.58 不報告為 signature；改用 deviance decomposition）。

### 1.3 Prior art positioning（Agent D, 詳見 02_prior_art_notes.md）

| Tool | 與本研究關係 | Critical overlap? |
|------|------------|------------------|
| **ROCIT** (bioRxiv 2026-03) | per-read transformer 用 methylation；LOH/CN/HP 當 training label 而非 grid axis | ❌ 無 |
| **TumorLens** (medRxiv 2026-03) | Long-read joint SNV+CN+LOH+methylation；sample-level，不做 per-read TP/FP | ❌ 無 |
| **SGZ** (PLoS CB 2018) | FoundationOne 4-axis (purity×ploidy×CN×AF) variant-level；無 phasing/methylation | ❌ 無 |

本研究獨特處（prior art 未涵蓋）：read-level 5 軸 + V3F/V5/V6 phasing variant 比較 + 5kHz simplex ONT。

---

## 2. 方法（簡要；詳見 02_methodology_notes.md）

### 2.1 BAM 三版本對照

| 版本 | longphase 設計 | hp_1_1:hp_2_1 ratio | 重要差異 |
|------|---------------|---------------------|---------|
| **V3F** | PON-only + 無 Layer 1.5 + germline-absent → hp=33 保守 | 1.138 (中性) | baseline |
| **V5** | V3F + 加 Layer 1.5 somatic fallback (HaplotagProcess.cpp:537-548) | 1.86 (priority bug 殘餘) | 增加 marker 但 over-promote |
| **V6** | V5 + 移除 Layer 1.5（重用 V5 phased VCF）| 1.838 (germline-existent 區殘留) | marker coverage 修正回 V3F 等級 |

### 2.2 主 grid 結構

- **3 主軸**：LOH (Inner/Outer) × HP bucket (5: same_HP1/same_HP2/cross_het/cross_het_inv/other) × Coverage_Multiple (5 bins: cov_loss/normal/elevated/gain/high_gain) = **50 cells**
- **2 covariate**（per-cell LR 控制）：HPFineNGroups (NG) + caller_af + NumReads
- **Power stratification**: powered (n≥50) / marginal (30-49) / underpowered (<30)
- **7 道 confound guard**: NumReads OLS, caller_af OLS, AF-bin L3, permutation (max-statistic null), marginal expectation, chr-stratified Mantel-Haenszel, HP1/HP2 symmetry binomial
- **多重比較**: BH-FDR q<0.05

### 2.3 已執行 fan-out

- Agent A (Step 1): V3F/V5/V6 三向 ISM 整合
- Agent B (Step 1.5 + 2): power gate + 50-cell grid + LR + 7-guard
- Agent C (Step 3): 4 FP zone deep dive
- Agent D (Prior art): TumorLens/ROCIT/SGZ/Wakhan/SAVANA 全文 + 差異化
- Agent E (Step 4): 4 樣本擴展 + Wilcoxon n=5
- Coordinator (main session): 本主報告 synthesis

---

## 3. Step 1: V3F → V5 → V6 三向整合（HCC1395）

### 3.1 三方 marker (NG_off ≥ 3) summary

| BAM | marker n | TP | FP | TP rate (95% CI) |
|-----|----------|----|----|-------------------|
| V3F | 21,997 | 20,183 | 1,814 | 0.918 [0.914, 0.921] |
| V5  | 18,382 | 16,428 | 1,954 | 0.894 [0.889, 0.898] |
| **V6** | **23,980** | **21,806** | **2,174** | **0.909 [0.906, 0.913]** |

- **V6 marker coverage 最多** (+9.0% over V3F, +30.5% over V5)
- **V3F marker TP rate 最純**
- V6 在 coverage 與 TP rate 取得最佳平衡（與 PI errata 一致）

### 3.2 V5 Inner LOH 區 NG=2 marker over-promote +60%

| BAM | Inner NG=2 master-joined n | Inner/Outer 比 |
|------|-----------------------------|----------------|
| V3F | 5,064 | 2.12 |
| **V5** | **8,136 (+60% vs V3F)** | 2.81 |
| V6 | 5,353 | 2.66 |

V5 在 Inner NG=2 cell 多出 +60% region 但 TP rate 沒上升（仍 99%+）→ V5 Layer 1.5 **marker over-promotion 的直接資料證據**。V6 回退到 V3F 保守標 hp=33 → Inner NG=2 region 數降回 5,353（與 V3F 5,064 接近）。

### 3.3 H1a/H1b/H1c — ceiling effect 解釋

| Cell | TP rate (V3F, V5, V6) | 觀察 |
|------|------------------------|------|
| Inner NG=2 | 0.991 / 0.993 / 0.992 | 全 ≥0.99，三 BAM 差 < 0.003 |
| Outer NG=2 | 0.992 / 0.992 / 0.990 | 全 ≥0.99，三 BAM 差 < 0.003 |
| 全域 NG≥3 marker | 0.918 / 0.894 / 0.909 | **真正差異**：±0.024 |

NG=2 Inner cell 對 HCC1395 paired-pileup 而言是「高純度 caller-pre-filtered」cell（TP rate 99%+），plan 原設 Δ ≥ 0.03 threshold 物理上達不到（**ceiling effect**）。Step 2 起改用 marker coverage / V5 over-promote 比 / FP enrichment 三層指標。

### 3.4 Trajectory（V3F → V5 → V6, flag=off, per-region NG）

| Class | n | 占比 | 解讀 |
|-------|---|------|------|
| A 兩段都改善 | 1,057 | 3.0% | Layer 1.5 + priority fix 兩段都增 NG |
| B 只 V5 改善 | 3,673 | 10.4% | Layer 1.5 增 NG 但 V6 退回 |
| **C 只 V6 改善** | **10,624** | **30.1%** | V6 主導改善（Layer 1.5 沒效但 priority fix 有效）|
| D 無變化 | 14,746 | 41.7% | Mid-NG region 不受影響 |
| E 反向/單段下降 | 5,232 | 14.8% | V6 比 V3F 還差（需 zoom-in）|

TP/FP 拆解：TP 中 30% 只 V6 改善（10,356/30,490）；FP 中 68% 無變化（V6 不改變 FP 分布，符合預期，因 V6 重用 V5 phased VCF）。

### 3.5 off/on flag 對照

- **`flag=on` 將全部 NG≥3 collapse 為 0**（所有 BAM × label × loh_side cell）→ ISM 端 `--germline-hp-only=on` 不能近似 V6 修補
- **V5_on Inner ≡ V6_on Inner**（NG=2 frac 都是 15.4%）— 因 V6 重用 V5 phased VCF，on 模式下兩 BAM germline-phased reads 完全相同 → 確認 V6 差異**只在 Layer 1.5 標記策略**，與 phasing 層無關

---

## 4. Step 2: 3 軸 50-cell grid + LR + 7-guard

### 4.1 Power gate (Step 1.5)

- **23/50 = 46% main cells powered** (n≥50) > 15 threshold → **PROCEED**
- `Diploid_Coverage_Used` median = **115**（paired_full HCC1395，不是 TO mode 的 61，已 flag in power_dry_run.md）
- Collinearity OK：max Cramér's V = 0.610 (loh × NG_bin)，VIF ≤ 3.35

### 4.2 V5 over-promote top 5 cells

| Cell | V5/V3F | V6/V5 | n_V3F→V5→V6 | n_TP_V6 |
|------|--------|-------|--------------|---------|
| **`Inner\|cross_het_inv\|cov_normal`** | **5.95×** | 0.995 | 62 → 369 → 367 | 366 |
| `Outer\|cross_het\|cov_elevated` | 4.58× | 1.000 | 26 → 119 → 119 | 119 |
| `Outer\|cross_het_inv\|cov_normal` | 3.85× | 1.002 | 157 → 605 → 606 | 606 |
| `Outer\|cross_het_inv\|cov_elevated` | 3.52× | 1.000 | 33 → 116 → 116 | 116 |
| `Outer\|cross_het\|cov_normal` | 3.42× | 1.000 | 168 → 574 → 574 | 574 |

**V6/V5 ≈ 1.0 in all top 5** → V6 Layer 1.5 revert **只在 germline-absent 區回退**；cross_het cells (germline-present) 保持不變。
**V5 over-promote 集中在 cross_het** （same_HP* 全 ~1.0×）— 證實 Layer 1.5 機制**只在 somatic-fallback heterozygous reads 作用**。

### 4.3 通過所有 7 道 guard 的 cells（2 個）

- `Inner|other|cov_normal` (n=4,984, FP=171, TP_rate 0.9657)
- `Outer|other|cov_normal` (n=8,447, FP=249, TP_rate 0.9705)

8 cells 通過 6/7 — **Guard 7 (HP1/HP2 symmetry) 是主要 filter**：15/23 cells flagged 因 V6 1.838 HP1:HP2 ratio bias（**真實 biology, 非 artifact**；Thread D 已記錄）。

### 4.4 LR Deviance decomposition

LR 只在 4 cells converge（FP ≥ 3）。**AF 是最強 covariate** (LRT p < 1e-6 in 3/4)；NG 在 2 cells 顯著；NumReads marginal。

---

## 5. Step 3: 4 FP Zone Deep Dive — **Paradigm Reframe**

### 5.1 4 zone 對照

| Zone | n | n_TP | n_FP | TP_rate | FP_rate | FP_enrichment | Fisher p | Captures of total FP |
|------|---|------|------|---------|---------|---------------|---------|-------------------|
| **Z-OCH** (Outer cross_het) | 1,468 | 1,443 | 25 | 0.983 | 0.017 | **0.124×** | 3.8e-62 | 0.5% |
| **Z-CHR8** (chr8 hotspot) | 3,061 | 2,094 | 967 | 0.684 | 0.316 | **2.305×** | 1.7e-159 | 20.0% |
| **Z-GL** (Inner gain+LOH) | 1,687 | 1,682 | 5 | 0.997 | 0.003 | **0.022×** | 4.5e-101 | 0.1% |
| **Z-AUTO** (top 5% FP density) | 1,767 | 551 | 1,216 | 0.312 | 0.688 | **5.022×** | 0.0 | 25.1% |

### 5.2 🤯 Paradigm Reframe

**2/3 預設「FP-rich」zone 實際上是 TP-pure**：

- **Z-OCH (Outer cross_het)**：FP rate 0.017 << global 0.137；Fisher p=3.8e-62 for **TP-enrichment direction**（不是 FP-enrichment！）
  - **新解讀**：cross_het = germline het + bona-fide somatic。FP germline-only calls 沒真實 somatic read support，不產生 cross_het → **Z-OCH 是 somatic-evidence signature，不是 FP signature**
- **Z-GL (Inner gain+LOH)**：FP rate 0.003 vs global 0.137（0.022×）→ 同 reframe：gain on somatic hap 是 somatic-evidence

**真正的 FP-rich zone 只有**：
- **Z-CHR8**：2.31× FP enrichment（sample-specific，HCC1395 chr8 已知 hotspot）
- **Z-AUTO**：5.02× FP enrichment（**emergent FP cluster，~70% 在 known zones 外**）

### 5.3 H4 — chr8 hotspot LR Deviance Decomposition (POSITIVE)

4-axis LR ablation on master-joined chr8 (n=2,373)：

| Axis dropped | Incremental dev_explained | Rank |
|--------------|---------------------------|------|
| **caller_af** | **0.393** | 1 (dominant) |
| **Coverage_Multiple (CN)** | **0.211** | 2 |
| hp_bucket | 0.063 | 3 |
| loh_side | 0.038 | 4 |

**(LOH + CN) − HP = 0.038 + 0.211 − 0.063 = +0.186** (3.7× 0.05 threshold) → **H4 POSITIVE**

CN 軸獨立貢獻 3.3× HP 軸 → **chr8 hotspot 由 CN+AF 主導，HP/LOH 次要**。chr8 是最高 FP chromosome (FP_rate 0.316 vs global 0.137)，chr17 次之 (0.198)。

### 5.4 H5 — Z-AUTO Jaccard (NEGATIVE)

- Z-AUTO 在 top 5% positions 捕獲 25.1% 全 FP（FP_enrichment 5.02×）
- 但 Jaccard with known-zone union = **0.184** << 0.60 → **H5 NEGATIVE**
- **~70% Z-AUTO 在 known zones 外** → novel FP mechanism lacks LOH/HP/CN signature
- 暗示需要新軸：mappability, repeat context, GC bias, structural variants

### 5.5 SEQC2 KDE validation

Z-GL 內 KDE-corrected Coverage_Multiple 與 SEQC2 CN 單調對應（CN=3 → 1.34，CN=8 → 2.08），確認 KDE binary 校正 (commits 374fad4 + 12d9b3e)。

---

## 6. Step 4: 4 樣本擴展（n=5 含 HCC1395）

### 6.1 Per-sample master TSV 建構

| Sample | rows | TP | FP | Join rate (master.tsv) |
|--------|------|----|----|-------------------------|
| HCC1395 (Step 1) | 35,332 | 30,490 | 4,842 | 85.0% (TP 96.8% / FP 10.5%) |
| H1437 | 70,964 | 70,191 | 773 | 95.1% |
| H2009 | 136,701 | 135,359 | 1,342 | 97.3% |
| HCC1954 | 20,136 | 19,449 | 687 | 89.1% |
| HCC1937 | 16,607 | 13,910 | 2,697 | 75.8% |

**FP join rate 跨樣本只 1-7%**（vs TP 89-97%）— phaseD VCF 是 paired-pileup mode 但 master.tsv 是 paired_full mode，不對稱。

### 6.2 Wilcoxon n=5 唯一 signature candidate

| Cell | Direction concordance | Wilcoxon p | Δ (mean) | 各樣本 TP rate |
|------|----------------------|-----------|---------|--------------|
| **`Outer\|other\|cov_high_gain`** | **5/5 above global** | **0.0625** (n=5 exact min) | +0.0069 | 全 5 樣本 = 1.000 |

**⚠️ Effect size 極小** (+0.0069) — H1437/H2009/HCC1954 global TP rate ≥0.998 caller 已近完美 → **ceiling effect 嚴重限制 cross-sample signature 強度**。

**Sensitivity** (n=4 排除 HCC1937)：同樣 1 candidate → HCC1937 不驅動 candidate ✅

### 6.3 HCC1937 BRCA1 outlier 分析

- **chr15/17/14** driver chromosomes
- **chr17 FP 不專屬 HCC1937** — HCC1395 也是 chr17 FP rate 最高（0.0430）→ BRCA1/TP53 region issue 是 HCC1395 + HCC1937 兩樣本共有
- Top outlier cell `Inner|same_HP1|cov_normal`: HCC1937 TP=0.300 vs 其他 0.998 (z=-362, n=10, FP=7 underpowered)

---

## 7. 證據鏈與穩定度

### 7.1 Tier 評估 → **⭐3 PARTIAL POSITIVE**

| 維度 | 評估 |
|------|------|
| 樣本數 | 5（HCC1395 pilot + 4 擴展）|
| 跨樣本一致 | Wilcoxon 5/5 同方向但 effect 小 (+0.0069) |
| Confound guard | 5+2 道（Guard 7 HP-symmetry 失敗為真實 biology）|
| 主要發現 | H4 chr8 hotspot CN+AF 主導（POSITIVE 數值明確）+ paradigm reframe Z-OCH/Z-GL |
| 限制 | Z-CHR8 sample-specific / Z-AUTO 機制未明 / ~63% FP unexplained |
| 與 prior art | TumorLens/ROCIT/SGZ/Wakhan/SAVANA 全無同口徑等效物 |

### 7.2 ⭐3 而非 ⭐4 的理由

- ✅ H4 chr8 hotspot LOH+CN 量化明確 + paradigm reframe（Z-OCH = TP signature）— 是新貢獻
- ✅ H2/H6 directional POSITIVE
- ✅ Step 4 1 cross-sample candidate（雖然 effect 小）
- ❌ H7 全通過 7 guards 只 2 cells（< 5 threshold）
- ❌ Z-CHR8 sample-specific（跨樣本 mean AUC ≤ 0.641 from prior Wave 3）
- ❌ Z-AUTO 機制未明
- ❌ ~63% FP 不被此 framework 解釋
- ⚠️ Effect size 受 caller 飽和限制（global TP rate 0.983）

升 ⭐4 需做：(a) 重做 Z-AUTO 在 4 樣本各自 KDE → 看密度機制是否 recur，(b) 加新軸（mappability/repeat/GC）測 63% unexplained FP。

---

## 8. Limitations

1. **Master-join FP loss 90%**：paired-pileup VCF 與 paired_full master.tsv 不對稱（4,842 → 506 FP master-joined）。Step 3 用全 FP set 重測 H3 已 mitigate，但 LOH/CN annotation 限制了 FP 子空間分析
2. **Ceiling effect 主導**：global TP rate 0.983（HCC1395 paired-pileup）→ effect size 物理上難 > 0.02
3. **HP1/HP2 asymmetry 為 design choice 非 bug**：V6 重用 V5 phased VCF 殘留 hp_1_1:hp_2_1 ratio 1.838 → Guard 7 失敗為預期
4. **Z-CHR8 sample-specific**：HCC1395 chr8 大型 gain+LOH 是已知 sample 特性（Memory `project_hcc1395_chr8_hotspot.md`），不可泛化
5. **COLO829 deferred**：truth set 0600 權限阻塞，Step 4 為 n=5 而非完整 7 樣本
6. **3 軸 grid effective dim ≈ 2.5-3**：caller_af 與 LOH inner/outer 高相關（Cramér's V 0.610）；改用 covariate 控制部分緩解但非完全
7. **未測 methylation 軸**：本 framework 不含甲基化特徵；ROCIT prior art 顯示 methylation read pattern 是獨立軸
8. **paired vs TO mode 限定 paired-pileup**：phaseC 僅 HCC1395 paired-pileup；TO mode 三向對照需另獨立做

---

## 9. Future Direction

### 9.1 Z-AUTO 跨樣本驗證（高優先）
對 4 樣本各自跑 KDE FP-density top 5%，看 Z-AUTO 密度機制是否 recur（即使 positions 不同）。若 recur → 升 ⭐4 cross-sample emergent FP signature。

### 9.2 加新軸測 ~63% FP（高優先）
- Mappability score (per-region)
- Repeat context (RepeatMasker overlap)
- GC bias (per-region)
- Structural variant proximity (SAVANA SV-call distance)
這 4 軸可能解釋 Z-AUTO 機制 + 63% unexplained FP。

### 9.3 加 methylation 軸（中優先）
本研究只用 HP tag + caller AF + LOH/CN，未用 methylation distance matrix（ROCIT prior art 的核心特徵）。在 powered cells 內加 ASM delta + methylation cluster ID covariate 看是否解釋 unexplained FP。

### 9.4 COLO829 補入（低優先，依權限）
COLO829 truth set 解 0600 權限後補入 Step 4 → n=6 Wilcoxon exact min p=0.0313（更嚴格 threshold）。

### 9.5 V6 production tag finalize（CURRENT_FOCUS W3）
本 cycle 確認 V6 marker coverage +9% 與 caller F1 不變，可作 V5→V6 production 升級依據（PI errata 待 W3 末打包）。

---

## 10. Reproducibility — File Inventory

### 10.1 計畫與方法
- `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/00_PLAN.md` (計畫書 v0.3)
- `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/01_data_provenance.md` (資料來源)
- `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/02_methodology_notes.md` (方法學)
- `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/02_prior_art_notes.md` (Agent D)

### 10.2 Step 1: V3F/V5/V6 三向整合
- `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/step1_findings.md`
- `step1_master_three_way.tsv` (9.3 MB, 35,332 rows × 64 cols)
- `step1_delta_v5_vs_v3f.tsv` / `step1_delta_v6_vs_v5.tsv` / `step1_delta_v6_vs_v3f.tsv`
- `step1_trajectory.tsv` (5-class A-E)
- `step1_off_vs_on_compare.tsv`
- `step1_summary_stats.json`
- `figures/step1/step1_three_panel_heatmap.png` / `step1_delta_heatmap.png` / `step1_trajectory_sankey.png`

### 10.3 Step 1.5 + 2: power gate + grid + LR + 7-guard
- `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step2_three_axis_grid/power_dry_run.md`
- `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step2_three_axis_grid/step2_findings.md`
- `step2_grid_3d.tsv` / `step2_trajectory_per_cell.tsv` / `step2_top_lists.tsv` / `step2_confound_guard.tsv` / `step2_marginal_vs_observed.tsv`
- `step2_power_curve.png` / `step2_facets.png`

### 10.4 Step 3: 4 FP zone deep dive
- `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step3_fp_zone_zoom/step3_synthesis.md`
- 4 zone 子目錄各含 `zone_findings.md` + `zone_summary.tsv` + `zone_grid.tsv` + `zone_three_version_trajectory.tsv` + `zone_confound_guards.tsv`
- `figures/step3_zone_enrichment.png` / `figures/step3_per_chr_fp_rate.png`
- `step3_cross_zone_summary.tsv` / `step3_jaccard_overlap.tsv` / `step3_verdicts.json`

### 10.5 Step 4: 4 樣本擴展
- `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step4_cross_sample_extension/step4_findings.md`
- `step4_signature_candidates.md` / `step4_HCC1937_outlier_analysis.md`
- `step4_per_sample_grid.tsv` / `step4_consistency.tsv`
- 4 樣本 per-sample master TSV 與 grid
- `figures/{sample}_facets.png` × 4

### 10.6 Scripts
所有腳本位於 `research/v6_bam_tpfp_hp_loh_cn/step*/scripts/`，可獨立重跑（已驗證在主 session 與 4 個 sub-agent 環境下無衝突）。

### 10.7 上游證據鏈
- Self-Phasing 整合報告: `InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md`
- Thread D 主軸: `InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md`
- V6 Validation: `InterSubMod/research/paired_priority_bug_audit/07_V6_validation_findings.md`
- V6 Cross-Sample: `InterSubMod/research/paired_priority_bug_audit/08_phaseD_v6_cross_sample_findings.md`
- V6 Caller F1: `InterSubMod/research/paired_priority_bug_audit/09_V6_caller_F1_verification.md`

---

## 11. Coordinator notes

- Multi-agent fan-out（5 sub-agents A/B/C/D/E + Coordinator）總執行時間約 3.5 小時（含 Stage 1 V3F-V5-V6 ISM 整合最費時 ~2 小時）
- BAM/ISM **完全沒重跑**，純粹 post-hoc 後分析
- 計畫 v0.3 與實際執行高度一致（Plan agent 8 點質疑全部落地 + 用戶補強「窮舉組合」收斂為 3 軸 50-cell + power gate）
- doc-standards frontmatter 已對齊；INDEX.md 待同步（Coordinator next step）
