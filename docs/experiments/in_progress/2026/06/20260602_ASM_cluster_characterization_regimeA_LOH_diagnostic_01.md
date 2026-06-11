<!--
build_date: 2026-06-02
agent: ASM cluster characterization (regime-A credible probe + LOH diagnostic)
status: in_progress
report_class: characterization (G1 read-level ASM clustering + G2 LOH diagnostic)
audience: PI / lab member / self (future trace)
sample: HCC1395 paired_full (single-sample; Tier-A ceiling, M9 cross-sample N/A)
task_type: B validation (whole-genome HP-axis, single-sample) — Tier-A characterization only
partial_flag: SINGLE-SAMPLE (HCC1395 only) — no cross-sample generalization claim
parent_research_dir: InterSubMod/research/tsg_promoter_asm_reviewer/
data_sources: research/tsg_promoter_asm_reviewer/genome_survey_v2/regimeA_credible_probe.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/regimeA_hardening.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/regimeA_residual_controls.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/loh_diagnostic_classifier.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/credible_loci_annotation.tsv
figures: figures/fig1_asm_landscape.png,figures/fig2_gate_evidence.png,figures/fig3_loh_triptych.png,figures/fig4_regression_to_extreme.png
verdict: G1 = credible-subset POSITIVE-but-bounded (15/23 survive full collider battery, above germline-het baseline Cliff δ=0.37) — Tier-A single-sample, refines (not overturns) 5/31; G2 = LOH apparent-two-haplotype 72% self-phasing artifact / 18% candidate subclone / 9% CN-regression (validates user blind-spot hypothesis)
last_verified: 2026-06-02
report_template: structured characterization + scientific-rigor §2-§7 evidence tiers
-->

# ASM 甲基分群 characterization：regime-A credible 子集 + LOH 成因診斷（HCC1395 單樣本）

> **scope note（誠實前提）**：本輪原計畫假設大量 greenfield，實際盤點發現 `research/tsg_promoter_asm_reviewer/scripts/15-35` + 20+ obs JSON 已完成約 70%。故本輪 = **consolidation + 補決定性 gate + 全新 LOH per-locus 診斷**，不是重新發現。用戶 2026-06-02 確認 scope = 嚴謹 characterization verdict、G1 framing = 先探 regime-A credible 子集、本輪限 HCC1395 單樣本明標。

---

## 0. TL;DR（verdict-first）

**G1 — 「跟 BRCA2 一樣清楚的 somatic 甲基分群位點」**：
在 somatic-controlled HP-axis、avoid regression-to-extreme 的 **credible regime（高 cov + 中段 baseline + nonLOH）** 子集中，**確有真 somatic 甲基分群訊號**：23/62 通過 length-placebo collider cut；其中 **15/23 撐過完整 artifact battery**（M4c CpG-context 96% + 強化 random-anchor M8 70%），且 regime-A 整體 blind-ARI **顯著高於 germline-het baseline**（median 0.229 vs 0.038，Cliff δ=0.37，p=2.3×10⁻⁴）。**這 refine（非推翻）5/31「ASM non-discriminative/coverage-modulated」結論** — 全基因組仍被 artifact 主導，但乾淨子集有真訊號。⚠ **嚴格誠實**：① 單樣本 Tier-A 天花板（M9 cross-sample 不可能 → 不可宣稱 Tier-B/generalize）；② 強化 collider gate **刷掉 SOX2/HOTTIP/SDHAF1 等好看的癌症基因**（它們是 genomic-context collider，非 somatic 特異）；③ BRCA2 本身在 regime-A 內但**未過 placebo cut**（borderline collider）。

**G2 — LOH unphase/HP3 盲點（用戶核心洞見）**：
LOH 區「表觀雙 haplotype」位點（longphase-S 不驗證 LOH 而產生）**經診斷分類，72% 為 self-phasing artifact（無真甲基 cluster）、18% candidate subclone、9% CN/regression** — **直接驗證用戶假設**：LOH 雙 haplotype 多數是 phasing artifact，但少數（18%）是 LOH 內真 subclone 甲基訊號。提供三成因 IGV-style read-level 圖各一 exemplar。

---

## 1. 兩個目標 → 方法對應

| 用戶目標 | 本輪做法 | 證據檔 |
|---|---|---|
| **G1** 找像 BRCA2 一樣清楚的 ASM 分群 + 量化「清楚」+ 依哪種 tag 一致 | regime-A credible 子集走 M0/M1/M2/M3/M5/M8/M4c 完整 gate（HP-axis）| `regimeA_credible_probe.json` + `regimeA_hardening.json` + `regimeA_residual_controls.json` |
| **G2** unphase/HP3/遠 read + LOH 盲點診斷 | LOH-regime 位點 CNV×甲基聚集**診斷分類**（artifact/subclone/CN），非強制歸類 | `loh_diagnostic_classifier.json` |

> **為何用 HP-axis 而非 ALLELE-axis**（回答 G1.2「依哪種 tag 一致」）：HP1 vs HP1-1 是**同一 germline haplotype 被 somatic allele 切分**，by construction 控制 germline allelic baseline；ALLELE-axis（ALT vs REF）= 兩條 germline haplotype 比較 = baseline allelic confounded（5/31 已證 TP 11.1% < het-null 15.2%）。**故 HP-axis 上超過 placebo 的分群才 somatic-attributable**。[L1, 對齊 ledger#55 + feedback_asm_allele_axis_baseline_confound]

---

## 2. G1 — regime-A credible 子集：有真 somatic 甲基分群（Tier-A 單樣本）

### 2.1 子集定義與 gate 結果（全 L1 實測）

regime-A = sig(p<0.05) HP-axis × n_paired_cpg≥100 × |germ_β−0.5|≤0.3 × nonLOH（obs23 定義）= **75 位點**；可評估 62（13 因 coverage 不足 drop）。

| Gate | 結果 | 判定 | 來源 |
|---|---|:---:|---|
| **M0** observed-only distance（禁 impute，drop poorly-connected） | matrix 無 NaN | — | `30_regimeA_credible_probe.py` (L1) |
| **M1** blind-ARI（盲分群恢復 somatic/germline split，PRIMARY） | median **0.229**，max 1.0 | — | probe.json (L1) |
| **M8** length-split placebo（collider proxy） | real >> placebo **p=5.07×10⁻⁸**（n=51 pair），median placebo ARI **0.023** | ✓ | probe.json (L1) |
| **M5** coverage-effect 簽名 | Spearman(ARI, log n_cpg) ρ=**−0.18 (p=0.16 NS)** → **無 coverage artifact** | ✓ | probe.json (L1) |
| **M2** vs germline-het null（baseline-allelic，n=60，cov-matched） | regime-A med ARI **0.229 vs 0.038**；MW **p=2.3×10⁻⁴**；Cliff δ **0.368**；median diff CI **[0.048, 0.337]**（下界>0）；pass-rate ARI≥0.30 **43.5% vs 15%** | ✓ | `regimeA_hardening.json` (L1) |
| **M3** rarefied silhouette（降採樣等組大小 B=200） | regime-A pass **0.308 vs het-null 0.080**（3.8×）→ 非 read-count artifact | ✓ | hardening.json (L1) |
| **M4c** CpG-context-removed（drop 變異±20bp CpG） | **22/23 (96%) 存活** → 非 CpG-destruction/sequence artifact | ✓ | `regimeA_residual_controls.json` (L1) |
| **M8-strong** random-non-somatic-anchor placebo | **16/23 (70%) 通過**；median anchor ARI 0.008；**7 個 collider** | ⚠ | residual_controls.json (L1) |
| **完整 battery 存活** | **15/23**（M4c=OK AND M8-strong=OK） | — | residual_controls.json (L1) |
| **M9** cross-sample sign-consistency | **未跑 = 單樣本不可能** → Tier-B 結構性不可達 | ✗(N/A) | — |

**funnel：75 regime-A → 62 evaluated → 23 length-placebo pass → 15 全 battery 存活**（見 Fig2 Panel B）。

### 2.2 關鍵誠實點（防 overclaim）

1. **單樣本 Tier-A 天花板**：M9 cross-sample（red-team A7 single-sample circularity = 本專案歷次 positive 全滅根因）需 ≥2 樣本 tagged BAM，目前只 HCC1395 → **不可宣稱 Tier-B / generalize**。[L1 框架 spec]
2. **強化 collider 刷掉好看的基因**：length-split placebo 太弱（within-locus），random-anchor M8 揭露 7 個是 genomic-context collider，**正包含 SOX2(anchor 0.20) / HOTTIP(0.115) / SDHAF1(0.169)** — 這些癌症相關 gene hit 的分群在隨機鄰近位點也出現 = 區域 read-geometry artifact，**不能當 somatic-ASM 發現**。[L1 residual_controls + annotation]
3. **BRCA2（chr13:32315128）本身**：在 regime-A 內，blind-ARI=0.79，但 placebo=0.132 > 0.10 → **未過 collider cut**（與既有 pertag「ALT⊂HP1-1 耦合」collider 疑慮一致）。BRCA2 是 reviewer 答辯 anchor、survey rank 25，但**不是最乾淨的 cluster**。[L1 probe.json line 25]
4. **15 個全 battery 存活者**多為 lncRNA/pseudogene/非經典位點（HERC6 / CIB2 / TBC1D16 / ENTPD8 / 多個 LOC）；**僅 13/37 promoter-proximal（|TSS|≤2kb）、18/37 在 CpG island/shore** → credible cluster 不全落在功能性 promoter。[L1 annotation TSV]

### 2.3 全 battery 存活的 15 位點（最近基因）

HERC6(chr4, ARI 1.0, shore) · LOC101927914(chr7, 0.71, island) · LINC00689(chr7, 0.55) · ENTPD8(chr9, 0.75, shore) · LOC105379446(chr9, 0.81) · LOC100128210(chr11, 0.79) · MIR9-3HG(chr15, 0.62) · **CIB2(chr15, 0.67, island, TSS+13)** · LOC124903622(chr16, 0.74) · **TBC1D16(chr17, ×3 loci, 0.37-0.42)** · ABCD1P1(chr2, ×2, 0.67) · RNU6-417P(chr7, 0.46)。[L1 cross of residual_controls + annotation]
*（SOX2/HOTTIP/SDHAF1/PTPRN2-AS1 等已被 collider/CpG-context gate 排除，不列為存活。）*

---

## 3. G2 — LOH 表觀雙 haplotype 成因診斷（驗證用戶盲點假設）

### 3.1 背景（用戶核心洞見，L1 工具碼確認）

longphase-S **不驗證 LOH**（無 LOH detection；germline-absent+somatic ALT → H3 / 升級 H1_1），不像 longphase-TO V6 用 heterozygosity-ratio 偵測 LOH 後給保守 `HP:i:33`。[L1 `SomaticHaplotagProcess.cpp:310-527` + `longphase-to/HaplotagProcess.cpp:515-526`]。故 LOH 區（物理單 haplotype）仍出現 HP1/HP1-1「表觀雙 haplotype」。⚠ **用戶提的「強制歸類 unphase REF→HP1、HP3→HP1-1」= longphase-TO V5 Layer 1.5 做過 → 證實造成 4.19:1 偏 HP1 + marker −19%，V6 已 revert**；故本輪採**診斷分類**（判成因）而非強制歸類。

### 3.2 診斷結果（L1）

LOH-regime sig HP-axis 位點 **2170**（cnLOH 816 / gainLOH 1343 / lossLOH 11）；分層抽樣評估 79：

| 成因 | n | % | 詮釋 |
|---|--:|--:|---|
| **self_phasing_artifact**（ARI<0.30 或 placebo collider） | **57** | **72%** | 表觀雙 haplotype **無真甲基 cluster** = 同一 haplotype 被變異切分的 phasing artifact |
| **candidate_subclone**（ARI≥0.30 & placebo<0.10 & 中baseline & 高cov） | **14** | **18%** | LOH 內**真甲基 cluster**（9 cnLOH + 5 gainLOH） |
| **CN_regression**（極端 baseline & 低 cov） | **7** | **9%** | regression-to-extreme（baseline 壓極端→假大 Δβ） |
| intermediate | 1 | 1% | — |

→ **直接驗證用戶假設**：LOH 雙 haplotype **多數（72%）是 phasing artifact**，但**18% 是真 subclone 訊號** — 不能一概而論。

### 3.3 三成因 exemplar（IGV-style read-level，Fig3）

| 成因 | exemplar | ARI | CN | n_cpg | germ_β | 説明 |
|---|---|--:|---|--:|--:|---|
| candidate subclone | **chr11:2146993 (IGF2)** | 0.892 | cnLOH | 112 | 0.59 | LOH 內真甲基分群 ⚠ **IGF2=imprinted gene，須謹慎詮釋**（雖 HP-axis 控制 imprinting） |
| self-phasing artifact | chr18:64001548 | −0.079 | lossLOH | 35 | 0.74 | germline/somatic 無分離 |
| CN / regression | chr8:134018510 | 0.523 | gainLOH | 28 | **0.96** | placebo=0.142 collider + 極端 baseline |

---

## 4. 圖表

| Fig | 檔案 | 內容 |
|---|---|---|
| Fig1 | `figures/fig1_asm_landscape.png` | 全基因組 HP-axis ASM landscape（Manhattan）；灰=全 4,694 sig（artifact-dominated）/ 綠=regime-A credible 75 / BRCA2 星標 rank 25 |
| **Fig2** | `figures/fig2_gate_evidence.png` | **headline**：A) regime-A vs het-null ARI violin（M2 δ=0.37）B) gate funnel 75→62→23→15 + scorecard（⚠ 強化 collider 刷掉 SOX2/HOTTIP/SDHAF1）C) LOH 3-class 分布 |
| Fig3 | `figures/fig3_loh_triptych.png` | LOH 三成因 read-level 甲基矩陣 triptych（candidate/artifact/CN 各一 exemplar） |
| Fig4 | `figures/fig4_regression_to_extreme.png` | coverage vs |Δβ| 散點 + LOH 著色 → 低 cov 放大 |Δβ|（regression-to-extreme 視覺證據） |

---

## 5. 方法完整性與未做項（誠實邊界）

**已跑 gate**：M0 observed-only · M1 blind-ARI · M2 germline-het null（n≥50 分層）· M3 rarefied · M5 coverage 簽名 · M8 length-placebo + **M8-strong random-anchor** · M4c CpG-context。
**驗尺（meta-validation）**：`24_cluster_eval_core.py` 已驗 PC1（simulated Δβ ARI>0.5）+ NC1（random-split ARI<0.15）= 尺基本 validated；in-battery 控制涵蓋 NC2(het-null)/NC3(coverage)/NC4(CpG-context)/NC5(collider) 主要失效模式。
**未做（remaining for full Tier-A / Tier-B）**：
- **M9 cross-sample**（BLOCKING for Tier-B）— 需 COLO829 等第 2 樣本 tagged BAM（MM/ML + HP:Z），目前不存在。
- 完整 PC1-3 / NC1-7 formal harness（本輪用既有 PC1/NC1 + in-battery 控制替代）。
- M4a mask-swap / M4b missing-permute / M4d basecaller-context（M4d 需 re-basecalling，工程量大）。
- regime-A 之外的 ALLELE-axis 位點（已知 confounded，不在 scope）。

---

## 6. Verdict 與 tier

| 主張 | tier | 證據等級 | 邊界 |
|---|:---:|:---:|---|
| credible regime（HP-axis）有真 somatic 甲基 cluster，高於 germline-het baseline | ⭐3 | L1（M1/M2/M3/M5/M8/M4c 實測） | **單樣本 Tier-A；M9 未驗 → 不可 generalize** |
| 全 battery 存活 15 位點為單樣本候選；SOX2/HOTTIP/SDHAF1 被 collider 排除 | ⭐3 | L1 | 需 cross-sample 才升 Tier-B |
| LOH 表觀雙 haplotype 72% self-phasing artifact / 18% candidate subclone | ⭐3 | L1 | 單樣本；CN class 為 categorical inference（非實測整數 CN，obs25 caveat）|
| 全基因組 ASM 仍 artifact-dominated（5/31 結論成立，本輪 refine 非推翻）| ⭐3 | L1 | — |

> **與 5/31 收斂的關係**：5/31「ASM real but non-directional / non-discriminative / coverage-modulated」在**全基因組 / ALLELE-axis / 全 regime** 尺度上仍成立；本輪 refine = **credible regime × HP-axis × 全 collider battery 下，存在 15 個高於 baseline 的真 somatic cluster 候選**（單樣本）。兩者不矛盾：之前看不到是因為被錯軸 + artifact regime 掩蓋。**filter-NEGATIVE（甲基當 FP filter, ⭐2 L4 DEAD）未被本輪觸碰** — 本輪是 characterization，非復活 filter。

---

## 7. Reproducibility

| 階段 | script | output |
|---|---|---|
| G1 regime-A probe（M0/M1/M8/M5） | `scripts/30_regimeA_credible_probe.py` | `genome_survey_v2/regimeA_credible_probe.json` |
| G1 硬化（M2 het-null + M3 rarefied） | `scripts/36_regimeA_hardening.py` | `regimeA_hardening.json` |
| G1 殘差控制（M4c + 強化 M8） | `scripts/37_regimeA_residual_controls.py` | `regimeA_residual_controls.json` |
| G2 LOH 診斷分類 | `scripts/38_loh_diagnostic_classifier.py` | `loh_diagnostic_classifier.json` |
| annotation | `scripts/39_annotate_credible_loci.py` | `credible_loci_annotation.tsv` |
| 圖 1/4 | `scripts/40_asm_landscape_figs.py` | `figures/fig1,fig4` |
| 圖 2 | `scripts/41_fig2_gate_evidence.py` | `figures/fig2` |
| 圖 3 | `scripts/42_fig3_loh_triptych.py` | `figures/fig3` |

BAM: `/big7_disk/.../longphase_s/HCC1395_tagged.bam`（longphase-S, HP:Z 字串 tag）。CNV: `/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed`。複用既有：`src/core/StructureTest.cpp`（PERMANOVA gold engine）、`24_cluster_eval_core.py`（M0/M1 base + PC1/NC1）、`cluster_quality_eval_framework.md`（M0-M9 spec）。

---

## 8. 後續（cross-sample 升級路徑）

1. **M9 cross-sample（最高優先，唯一通往 Tier-B）**：COLO829 / 其他 paired tagged BAM 跑同 15-locus battery，檢 sign-concordance。
2. 完整 PC1-3 / NC1-7 formal meta-validation harness。
3. M4a/M4b/M4d 補完整 Tier-A battery。
4. 15 存活位點的功能跟進（motif / expression）— 但須先 cross-sample 確認非單樣本 artifact。
