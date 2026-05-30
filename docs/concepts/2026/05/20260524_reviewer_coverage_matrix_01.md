---
date: 2026-05-24
topic: LongPhase-S reviewer comments × ISM deliverable A coverage matrix
source_files:
  - InterSubMod/../../big8_disk/liaoyoyo2001/knowledge/raw/paper-notes/longphase-s-reviewer-response-2026-05-23/report.md (291 lines)
  - InterSubMod/docs/reports/in_progress/2026/05/20260523_HKU_collab_tag_proposal_01.standalone.html (1167 lines)
  - InterSubMod/.claude/memory/project_loh_subclone_af_methylation_positive.md (31d old)
  - InterSubMod/.claude/memory/project_loh_constrained_phasing_discovery.md (31d old)
  - InterSubMod/.claude/memory/project_hpfinengroups_subclone_marker.md (12d old)
  - InterSubMod/.claude/memory/project_phase2_cycle1_global_fp_filter.md (3d old)
  - InterSubMod/.claude/memory/project_clairs_to_verdict_pilot.md (34d old)
framework: "Coverage Matrix"
purpose: "5/24 HKU ClairS-TO Luo Lab handoff alignment — map our ISM deliverable A evidence onto LongPhase-S reviewer comments to anticipate which reviewer concerns we can credibly address (vs honestly skip)"
audience: PI + HKU-BAL ClairS team
git_head: 369e8f9 (refactor/phase1-safety)
rigor_tier: Tier 3 (≥500 字, cross-file mapping, evidence-tier annotated per `/scientific-rigor` §2)
---

# LongPhase-S Reviewer × ISM Deliverable A Coverage Matrix

> ✅ **T10 全 24 chr 完成（2026-05-24 18:00 update）— HP3 雙峰發現升級**：原 R3-M5 row T1 結論「HP3 是 somatic-evidence bucket」修正為 bimodal — Type A (20/23 chr) HP3 frac 0.22-0.38% / TP rate mean 90.4% = somatic-evidence (6,900× enrichment)；Type B (chr6 HLA / chr16 segdup / chrX) HP3 frac 8.6-21.6% / TP rate 0-3.5% = phasing-failure fallback。Pooled HP3 1.76% / TP 14.3%。R3-M5 row 嚴重度升 ⭐⭐⭐⭐⭐（concept-leading）。完整: `InterSubMod/research/hku_collaboration/findings_5_24/T10_HP3_24chr_findings.md` + T11 4-layer figure `InterSubMod/research/hku_collaboration/figures/T11_4layer_seetree_and_forest.png`。方法論修正 per memory `feedback_observation_scope_default_comprehensive.md`。

## §1 TL;DR (v6 audit-corrected)

我們對 reviewer 23 條評論（R1×7 + R2×4 + R3×12，含 1 隱含 R3-Curiosity-Implicit）可形成的對應如下：**5 條強對應**（ISM 直接拿出 L1-L2 evidence figure/數據；可主動引用協助 LongPhase-S rebuttal）、**6 條部分對應**（有相關 evidence 但屬 L2-L3 觀察 / 不同 sample / 機制需 reframe，需加 caveat 才能引用）、**12 條應誠實 skip**（屬 LongPhase-S 本體 architecture / 純排版 / 純 caller mechanism / out-of-scope）。**Total: 5 + 6 + 12 = 23 entries**（v5 errata：原 TL;DR 寫 22 為計數錯誤，§2 table 與 §5 skip table 已對齊 23）。

強對應 5 條依強度排：**R1-Curiosity（TSG×methylation）** > **R3-M6（multi-subclone）** > **R1-M1（subclonal/region GHIR）** > **R3-M5（HP3 fraction）** > **R3-M1（rescued/removed variant region/class）**。其中 R1-Curiosity 是 reviewer 自行 cue 我們的延伸研究方向（Discussion 已預留 hook），對 5/24 handoff 是最自然的對話切入點。

最該 skip 的是 ASCAT/Battenberg/PURPLE benchmark（R2-#4）、Fig 排版 typos（R3 Editorial）、cellular vs DNA purity 概念辯論（R2-#1）、SNV-vs-indel mechanism（R3-M3）、callset overlap Venn（R3-Sp2）— 這些屬 LongPhase-S architecture / 排版 / 純 caller mechanism，ISM 介入只會稀釋論述。

## §2 完整對照表（23 entries 全列）

> 註：reviewer comments source 文件共 sections 1.1-1.7（R1, 7 items）+ 2.1-2.4（R2, 4 items）+ 3.1-3.12（R3, 12 items）= **23 enumerated**。表中 R3-Editorial（合併 3.11 typos + 3.12 figure detail）為 1 條 → **22 distinct issues**。「17 條」如指剝離掉所有純編輯類 + figure 排版 minor items，則為 R1×5 + R2×4 + R3×8 substantive = 17（在「嚴重度 ≥⭐2」欄一致）。

| ID | 評論摘要 (≤60 字) | 嚴重度 | 我們狀態 | Evidence 來源 (path:line / memory) | Tier [F/O/I/U] | 我們報告對應位置 | 回應建議 (≤80 字) |
|---|---|---|---|---|---|---|---|
| **R1-M1** | subclonal / region-level GHIR：實體腫瘤有 subclone+CNV+LOH，purity 應視為 *effective fraction*；region GHIR 可揭示 subclone | ⭐⭐⭐⭐ | ✅ 強對應 | `project_loh_constrained_phasing_discovery.md` (6/6 TO samples Inner same-hap 93-99%, TP gap +0.37) + HKU report §4.7 A8 ideogram (chr8 96% LOH) | [O] cross-sample observation | §4.5 A6 / §4.7 A8 / §1.4 A2 | 引用 HCC1395 chr8 99.1% LOH × HP2-1 dominance + Inner/Outer NG=2 split 為 "region-level GHIR-equivalent subclone signature" precedent |
| **R1-M2** | read-level benchmark included-read fraction 未報告 | ⭐⭐⭐ | ⚠ 部分對應 | HKU report §1.4 A2_2 stat-grid (`pair common reads 70.5 within ±16 bp, 16 over 500 bp`) + A2_1 (85.9% PS ≥2 TP, n=759) | [O] HCC1395 single-sample | §1.4 A2_1 / A2_2 | 提我們 ISM 也維持「read carrying ≥2 anchor」coverage 統計慣例；不直接幫 LongPhase-S 算 fraction 但 demo 統計範式 |
| **R1-m1** | haplotype-specific CNV × GHIR 互動分析缺失 | ⭐⭐ | ⚠ 部分對應 | HKU report §4.5 A6 6-zone stratification (Kruskal-Wallis H=234.11, p=1.41e-48) + §4.7 A8 (chr20 94.5% CNV / 0% LOH HP1:HP2 對稱) | [O] HCC1395 ONT only | §4.5 A6 / §4.7 A8 | 引 A6 結論 "subclone proxy 不能由 LOH×CNV stratify" + A8 "LOH ≠ phasing failure" — 對齊但 caveat：cell line + single platform |
| **R1-m2** | tumor / normal coverage range 沒寫進正文 | ⭐ | ❌ skip | — | — | — | 排版細節；LongPhase-S 自己改即可 |
| **R1-m3** | somatic haplotype block size 分布未報 | ⭐⭐⭐ | ⚠ 部分對應 | HKU report §4.6 A7 (Global N50 943 kb, max 18.3 Mb chr1, 6,861 blocks) | [F] direct measurement from PS tags in HCC1395 BAM | §4.6 A7 | 提我們 PS block N50 量化方法（per-chr ranking, chr22=1.625 Mb max autosome, chrX=101 kb min）作 LongPhase-S 自報範本 |
| **R1-m4** | Fig 6 legend × page number 重疊 | ⭐ | ❌ skip | — | — | — | 排版；無關 ISM |
| **R1-Curiosity** | TSG promoter hypomethylation × somatic haplotype（非強制；Discussion 預留 InterSubMod hook） | ⭐⭐⭐⭐⭐ | ✅ 強對應 (signature) | `project_loh_subclone_af_methylation_positive.md` (ΔNG=+0.787 7/7 POSITIVE, p<10^-65) + `project_hpfinengroups_subclone_marker.md` (chr19 Phase B 94.7% TP marker rate) + HKU report §2 read×CpG matrix design + `InterSubMod/src/core/MatrixBuilder.cpp` + `DistanceMatrix.cpp` | [O] 7-sample paired + cross-binary; [I] mechanism reframe LOH-constrained phasing | §2.1-§2.5 / F5 ISM heatmap | **這就是我們的 hero deliverable** — 提 "MethylSomaticAnalysis / InterSubMod follow-up study on TSG ASM × HP1-1/HP2-1 sub-clone branches"；可挑 TP53/CDKN2A/RB1 之一做 single-region IGV demo |
| **R2-#1** | cellular purity vs DNA fraction (ploidy/WGD 混淆) | ⭐⭐⭐ | ❌ skip | — (這是 LongPhase-S 的 purity model 論述問題，與 ISM HP tag downstream 無關) | — | — | 不主動回應；屬 LongPhase-S in silico mixing 設計 + per-cell-line label，ISM 介入無 leverage |
| **R2-#2** | 訓練/測試 overlap → overfitting (事實錯誤) | ⭐⭐ | ❌ skip | — | — | — | LongPhase-S 自行引 Methods 4.3.2 + leave-k-out CV R²=0.96 駁回；ISM 無 standing |
| **R2-#3** | 訓練樣本過 extreme，需更多 low-rearrangement cancer genome | ⭐⭐⭐⭐ | ⚠ 部分對應 (warning) | `project_phase2_cycle1_global_fp_filter.md` (HCC1395 in-distribution 5-fold OOF +0.02236 / HCC1395 LOSO held-out **−0.00012** = 100% effect circularity / 5-sample LOSO mean **−0.00004** / HCC1954 transfer **−0.377** catastrophe) + `project_clairs_to_verdict_pilot.md` (purity=0.40 subsample F1 gain=0 與高純度產線外推差) | [F] **negative**: ISM filter 自己也 fail cross-sample LOSO | §4 Q4 (Phase 2 status note) | 誠實提我們 ISM filter 也遇 cross-sample failure → "cell-line dataset 外推到 patient sample 是 common pitfall，跨 5 樣本 LOSO mean ≈0 DIRECTION_NEGATIVE Wilcoxon p=0.125" — 用 negative case study 而非 over-claim |
| **R2-#4** | ASCAT 跑法不透明 + Battenberg/PURPLE 應跑 | ⭐⭐⭐ | ❌ skip | — | — | — | 純 LongPhase-S benchmark methodology；ISM 無對應 baseline |
| **R3-M1** | recalibration 後哪類 variant 救回/砍掉（region/class/GC/VAF bin） | ⭐⭐⭐⭐ | ✅ 強對應 | HKU report §3.2 F7 chr2:18,086,020 14-site stratification (4-Layer Evidence: SEQC2 truth + LOH context + HP2-1 dominance + caller blind spot) + `project_loh_constrained_phasing_discovery.md` Inner vs Outer NG=2 TP rate (0.93 vs 0.55) | [O] HCC1395 single-region pilot + 6-sample TO cross-sample | §3.2 F7 + §1.4 A2_2 | 引 F7 "LOH × HP2-1 集中 region" 為 rescued variant class 範例 + Inner same-hap NG=2 TP enrichment 為「LongPhase-S HP1-1/HP2-1 救回的 variant 屬何類」的機制詮釋 |
| **R3-M2** | PacBio HiFi 跨平台未測 | ⭐⭐⭐ | ❌ skip | — (我們手上 HCC1395 / 7 sample 全 ONT) | — | — | 平台覆蓋是 LongPhase-S 自己的義務；ISM HCC1395 = 5kHz ONT 一致 |
| **R3-M3** | SNV vs indel benefit 差異未解釋 | ⭐⭐ | ❌ skip | — (我們 ISM 不處理 indel) | — | — | 純 caller mechanism；不主動回應 |
| **R3-M4** | filter threshold 5 purity bin 數值未列 | ⭐⭐ | ❌ skip | — | — | — | LongPhase-S 自身 Methods 4.4.5；無 ISM 連結 |
| **R3-M5** | HP3 佔比 + HP3 中 TP rate 未報 | ⭐⭐⭐⭐⭐ | ✅ 強對應 (concept-leading) | **T10 (5/24 全 24 chr) HP3 bimodal**: Type A 20/23 chr HP3 frac 0.22-0.38% / TP rate 78.6-95.4% mean 90.4% (= somatic-evidence bucket, 6,900× random enrichment)；Type B 3 chr (chr6 HLA 9.41%/1.37% + chr16 segdup 8.58%/3.53% + chrX 21.60%/0% truth=0) = phasing-failure fallback。Pooled (24 chr) HP3 1.76% / TP 14.30% (outliers 主導 85.4%)。HKU report §1.2 HP 5-階層 table + §4.7 A8 per-chr HP summary | **[F]** direct measurement T10 (24-chr scan, sanity check chr1/chr8/chr19 對齊 T1 baseline) | §2.5 + §1.2 / §4.7 A8 + T11 4-layer figure | 提我們 T10 bimodal finding + T11 4-layer integrated figure (aggregate scatter / canonical chr1 / extreme chr6+chr16+chrX outlier / well-explained F7 chr2:18,086,020 14-site) 為 LongPhase-S HP3 prototype panel template。雙重定義避免單一過度宣稱：well-phased 區 somatic-enriched / HLA-segdup-sex 區 fallback |
| **R3-M6** | multi-somatic-haplotype（HP1-2 / HP1-3）未涵蓋；IGV 範例 | ⭐⭐⭐⭐⭐ | ✅ 強對應 (concept-leading) | HKU report §2.3 HP 細分樹 linear depth notation (HP1-1-1 / HP1-1-2 / HP1-2-1) + F4v2 4 scenarios (A linear / B branching / C branching+linear / D HP1+HP2 並行) + `project_hpfinengroups_subclone_marker.md` ⭐3 pipeline-dependent marker | [I] design proposal; [O] 4-binary cross-sample marker rate 0.817-0.992 (5/7 ≥0.85) | §2.3 F4v2 + §1.2 table | **這就是 deliverable A 的 core proposal** — 我們直接給出 LongPhase-S HP1-2 / HP1-1-1 linear depth schema + R/A/X 9-state matrix + X 缺失 ISM 補完規則 (RX→methyl cluster→vote)；F7 chr2:18,086,020 即可作為 multi-subclone IGV 範例（雖目前 HP2-1 dominant 不展示 multi-branch，但 frame ready） |
| **R3-M7** | 低 VAF precision/recall + detection 下限 | ⭐⭐⭐ | ⚠ 部分對應 (mechanism) | `project_hpfinengroups_subclone_marker.md` (NG=4 + AF<0.4 + NR≥80 NonLOH TP rate 0.928, AF<0.2 TP 0.937; HCC1954 AF≥0.4 FP enrichment fraction 0.70) + `project_phase2_cycle1_global_fp_filter.md` (caller_af coef +3.44 dominant) | [O] 7-sample HCC1395 + HCC1954 stratified | §4 Q4 (caller_af 5th-rank statement) | 提 AF<0.4 是我們 ISM filter 也驗證的 "low-VAF region somatic-rich" rule + HCC1954 AF≥0.4 FP catastrophic enrichment 可作 LongPhase-S 低 VAF panel 對照 (不同 sample 換 frame) |
| **R3-Sp1** | ClairS 比 DeepSomatic 受益更多原因 | ⭐⭐ | ❌ skip | — | — | — | 屬 caller training data + neural arch 差異；ISM 無 leverage |
| **R3-Sp2** | DeepSomatic / ClairS callset overlap Venn + 為何選 ClairS+LPS | ⭐⭐ | ❌ skip | — | — | — | 純 caller comparison；不主動回應 |
| **R3-Sp3** | germline phasing tool（WhatsHap/LongPhase）head-to-head diploid vs LOH 區段 | ⭐⭐⭐⭐ | ⚠ 部分對應 | `project_loh_constrained_phasing_discovery.md` (LongPhase-S `--somaticMode` 在 LOH Inner 區 NG=2 = 93-99% same-hap, 對比 germline phasing 預期 50/50 cross-het) + HKU report §1.4 reframe (germline het 不變但 HP imbalance shifts) | [O] 6-sample TO cross-sample mechanism | §1.4 + Q4 reframe | 提 "LongPhase-S `--somaticMode` vs mainline germline-only phasing" 機制差 — 但我們未跑 LongPhase 原版對照，僅有 V5/V6 binary 比較 (`paired_priority_bug_audit/08_phaseD`)，純 same-binary 不同 flag |
| **R3-Sp4 + R3-Fig** | typos (haplotye, hypotehsize, Uknown) + Fig 1B/2A/6B legend illegible + Fig 6A "Uknown" | ⭐ | ❌ skip | — | — | — | 純編輯類；無 ISM 介入空間 |
| **R3 Curiosity-Implicit** | (隱含) Purity scale 方向 / "Called Somatic Variants" 池 confusion | ⭐ | ❌ skip | — | — | — | Figure 排版；skip |

## §3 強對應 5 條逐條展開

### 3.1 R1-Curiosity — TSG promoter hypomethylation × somatic haplotype (我們最強對應)

**Reviewer 原文（翻譯）**：「Discussion 已提到 'integration of orthogonal signals such as methylation phasing may enhance...'，這剛好對應你們知識庫裡 MethylSomaticAnalysis / InterSubMod 的延伸主題。可以在 rebuttal 提一句 'we plan to demonstrate this in a follow-up study with MethylSomaticAnalysis / InterSubMod'，或挑一個熟知的 TSG（TP53, CDKN2A, RB1）做 IGV 範例放進 Supp。」

**我們的 evidence**：
- **Paired AF × NGroups POSITIVE**（`project_loh_subclone_af_methylation_positive.md`）：ΔNG=+0.787 跨 7 樣本 paired mode（p<10^-65），median |r|=0.755，6/7 segment ρ POSITIVE — methylation 多樣性與 LOH intermediate AF 強正相關。[O] tier。
- **HPFineNGroups chr19 conditional POSITIVE**（`project_hpfinengroups_subclone_marker.md`）：V6 chr19 Phase B cross-tab 768 regions，flag=off NG≥3 TP rate **94.7%**，flag=on NG_on=2 **91.5%**，最強 cell NG_off=5→NG_on=2 (122 regions) TP rate **99.2%**。[O] tier with caveat（chr19 only = 2.16% genome, BAM = V3F binary 非 V5）。
- **5/7 樣本 V6 cross-sample 驗證**（同 memory）：marker rate 0.992 / 0.993 / 0.954 / 0.817（4 樣本，COLO829 deferred），3/4 ≥0.85 gate。
- **ISM C++ 模組已實作**：`InterSubMod/src/core/MatrixBuilder.cpp:17-39`（read × CpG matrix builder, MM/ML tag 解析）+ `InterSubMod/src/core/DistanceMatrix.cpp:310-345`（6 距離度量 source-verified: **NHD (Hamming) / L1 (Manhattan) / L2 (Euclidean) / CORR (Pearson) / JACCARD / BERNOULLI**；**無 Cosine**，v5 errata 修正）— 可實機跑 single TSG region demo。

**在 HKU 對話中如何引用**：
> "對 reviewer 1 curiosity 的 TSG×methylation 設想，我們手上有 7-sample paired mode 的 ΔNG=+0.787 訊號 + V6 chr19 phase-resolved NG=5→2 cell TP rate 99.2% 的 marker — 屬 [O] 觀察級別。我們可在 LongPhase-S Supp 補一個 TP53 / CDKN2A / RB1 IGV 圖：(1) ClairS-TO ssrs tp.vcf 標 PASS variant；(2) LongPhase-S HP1-1/HP2-1 分支 read 顏色；(3) ISM read×CpG heatmap 顯示 5mCG 機率聚類。但須註明：(a) 此為 follow-up scope；(b) 不是 LongPhase-S 本身 endpoint metric；(c) HCC1395 / 7-sample 是 cell line。"

**Caveat**：別把 ΔNG=+0.787 直接稱「ASM rescue gain」— TO 層機制已 reframe 為 LOH-constrained phasing signature（`project_loh_constrained_phasing_discovery.md`），Paired 層加註「需獨立 phasing-vs-methylation 驗證」。

---

### 3.2 R3-M6 — Multi-somatic-haplotype（HP1-2 / HP1-3）+ IGV 範例

**Reviewer 原文（翻譯）**：「cancer cell line 應該會出現多個 somatic subclone（HP1-1, HP1-2, …），論文完全沒提。請說明這種情境下表現如何，或若被當 artifact 砍掉，請用一張 IGV（類似 Fig.1 LOH 例子）展示。」

**我們的 evidence**：
- **HKU 報告 §2.3 已給 linear depth notation**：`HP1-1`（HP1 第一代 child，+1 somatic event α）/ `HP1-1-1`（linear extension, +2 events）/ `HP1-2`（HP1 第一代 sibling，不同 mutation event）/ `HP1-1-2`（HP1-1 第二代不同 event sibling）。[I] design proposal。
- **F4v2 4 scenarios 圖**：A linear / B branching 1st-gen / C branching+linear (AA ambiguous 紫色) / D HP1+HP2 並行（germline het root）。
- **R/A/X 9-state matrix**（HKU 報告 §2.1）+ X 缺失 ISM 補完規則（§2.2 RX→同 PS read→methyl cluster→鄰近 SNV vote→{RR, RA}）。
- **F7 chr2:18,086,020 14-site stratification**（HKU 報告 §3.2）— 雖然目前 HP2-1 dominant（71%）僅展示 single sub-branch，但 frame 已 ready：14 個 anchor 位點 + LOH context + HP simbalance 證據。

**在 HKU 對話中如何引用**：
> "對 reviewer 3 Major #6 multi-subclone limitation，這正是我們 deliverable A 的 core proposal — 我們提議把 LongPhase-S 既有 5 階層 HP tag 擴展為 linear depth schema (HP:Z:1-1-1 / 1-1-2 / 1-2-1)，配合 ISM read×CpG methyl cluster 補完 X 缺失 — encode 多 subclone 演化樹進 BAM aux field。F4v2 圖 + R/A/X matrix 為 design proposal [I]；F7 chr2:18,086,020 為 single-region precedent。建議 LongPhase-S Methods/Discussion 引入這個 schema 作 'future extension to multi-subclone modeling'。"

**Caveat**：linear depth schema 是 design proposal [I]，未跑 multi-subclone benchmark；F7 region 目前只展示 1 個 sub-branch（HP2-1 dominance），未實證 HP1-2 / HP1-3 多分支 case。

---

### 3.3 R1-M1 — Subclonal / region-level GHIR

**Reviewer 原文（翻譯）**：「LongPhase-S 目前回報的是 genome-wide 單一 tumor purity；但實體 primary tumor 常有 subclone、CNV、LOH 造成 region-specific GHIR 偏移。請補述：(a) 多 subclone 情境下這個 purity 該被視為 *effective / averaged tumor cell fraction*；(b) 局部 GHIR 偏離 genome-wide 分布是否可進一步透露 subclonal 結構。」

**我們的 evidence**：
- **`project_loh_constrained_phasing_discovery.md`**：TO 模式跨 6 樣本 obs18 — Inner × NG=2 組成 93-99% same-hap（HP1+HP1-1 或 HP2+HP2-1）；Outer × NG=2 cross-het 占 0.1-97.5% median 44%；TP rate gap (Inner same_HP1 − Outer cross_het) median +0.37, 6/6 樣本一致正向。HCC1937 +0.52 / HCC1395 +0.46。[O] tier。
- **HKU 報告 §4.7 A8 24-panel chr ideogram**：chr8 96% LOH / chr17 93% / chr11 85% (Tier 1) vs chr16 0% / chr20 0% / chrX 0% (Tier 4 reference)。chr20 = 94.5% CNV coverage 但 0% LOH → CNV gain 不必伴 LOH，HP1/HP2 仍對稱 (28% / 27%)。
- **HKU 報告 §4.5 A6 6-zone stratification**：Kruskal-Wallis H=234.11 p=1.41e-48；方向與預期相反（LOH inside subclone ratio 0.07-0.09 反而最低，Z6 0.435 anomaly 經審查為 chrX X-inactivation + chr6/chr16 centromere artifact）。

**在 HKU 對話中如何引用**：
> "對 reviewer 1 Major #1 region-level GHIR，我們有 6-sample TO Inner same-hap NG=2 = 93-99% / TP rate +0.37 cross-sample 證據（[O] 觀察）— LOH region 物理上只剩單 haplotype，somatic SNV 必然產生 same-hap ref/alt 子族（=NG=2 same-hap）；非 LOH 區的 NG=2 多為 cross-hap canonical het phasing。建議 LongPhase-S Discussion 引入 'Inner same-hap NG=2 signature' 作為 region-level subclone evidence 範例。我們的 A8 ideogram 可作 per-chr LOH × HP composition 對照 supp figure 範本。"

**Caveat**：A6 顯示「subclone proxy 不能由 LOH×CNV stratify」，建議避免直接宣稱「ISM 給出 region GHIR」— 改提 "我們 quantify HP composition shifts at chr-resolution but mechanism is phasing not methylation per se"。

---

### 3.4 R3-M5 — HP3 佔比 + HP3 中 TP rate

**Reviewer 原文（翻譯）**：「論文沒講 HP3 佔總 variant calls 的比例，也沒講 HP3 之中 TP 的比例。」

**我們的 evidence**：
- **HKU 報告 §1.2 HP 5 階層 table**：HCC1395 V6 全基因組 autosome 平均 HP3 < 1% per autosome（[F] from V6 tagged BAM direct count）。
- **HKU 報告 §4.7 A8 per-chr HP summary**：chr6 outlier HP3+HP1-1+HP2-1 合計 ~9%；chr16 no_HP 52.9% (autosome 最高)；chrX HP3 12.0% / no_HP 61.1% (single-X 簽名)。
- **HKU 報告 §4.6 A7 N50**：chr6 N50 734 kb / chr16 226 kb / chrX 101 kb — 與 A8 HP3+no_HP 比例 cross-validate（同一軸向兩種量化）。

**在 HKU 對話中如何引用**：
> "對 reviewer 3 Major #5 HP3 fraction，我們的 24-panel ideogram (HKU 報告 §4.7 A8) per-chr HP1/HP2/HP1-1/HP2-1/HP3/no_HP normalized stack 即是 LongPhase-S HP3 panel 的範本（autosome 平均 <1%, outlier chr6 chr16 chrX 標出）。建議 LongPhase-S Supp Table 給每 dataset 的 (HP1-1 / HP2-1 / HP3 / no_HP read count, fraction, F1 in subset)。我們 HCC1395 V6 BAM 數據可分享作 prototype。"

**Caveat**：HP3 中 TP rate 我們未直接算 — 我們的 F1 統計是 caller-level（ClairS-TO ssrs）而非 HP3-stratified。

---

### 3.5 R3-M1 — 救回 / 砍掉的 variant 屬何 class / region

**Reviewer 原文（翻譯）**：「希望看到 LongPhase-S recalibration 後，哪類 variant 變化最大（重複區、TR、SD、GC bias）？被『救回』的 variant 為什麼一開始被原 caller 漏掉？以及 callset 大小前 / 後對照。」

**我們的 evidence**：
- **HKU 報告 §3.2 F7 chr2:18,086,020 14-site stratification**：4-Layer Evidence — Layer 1 SEQC2 truth (125 TP/Mb vs chr2 avg 10 TP/Mb), Layer 2 LOH context (14/14 in LOH), Layer 3 HP2-1 dominance (58-76% avg 67%), Layer 4 caller blind spot (nearest in_loh anchor 4,332 bp >> 16 bp window)。
- **`project_loh_constrained_phasing_discovery.md`**：Inner same-hap NG=2 TP rate 0.93 vs Outer cross-het 0.55 — gap +0.37 跨 6 樣本。即 "LOH-driven sub-clone TP candidate" 可從 "true FP" 區分。
- **`project_hpfinengroups_subclone_marker.md`**：HCC1954 / HCC1937 AF≥0.4 FP enrichment fraction ~70%（vs AF<0.4 normal）+ all 5 in-scope 樣本 clustering ratio <1（非 spatial 而是染色體層級 CNV-driven）。

**在 HKU 對話中如何引用**：
> "對 reviewer 3 Major #1 rescued / removed variant class，我們的 F7 chr2:18,086,020 32 kb region 是 single-region pilot：LOH × HP2-1 集中區，14 anchor sites avg HP2-1 67%。這代表 LongPhase-S '救回的 TP candidate' 屬 'LOH-driven sub-clone HP-imbalance signature' 一類。HCC1954 AF≥0.4 FP enrichment 0.70（vs AF<0.4 normal）為 '砍掉的 FP' 對照範例。建議 LongPhase-S Supp 加 region/class stratification table（LOH yes/no × HP imbalance yes/no × AF bin × variant class）。"

**Caveat**：我們的 region stratification 是 per-sample post-hoc，非 LongPhase-S 端 native；只能提供 framework 而非 LongPhase-S 直接結果。

## §4 部分對應 + 我方 caveat 列表

| ID | 對應強度 | 我方 caveat |
|---|---|---|
| **R1-M2** (included-read fraction) | 部分 (provide stats范例) | A2_2 fraction 是 HCC1395 single-sample [O]，非 LongPhase-S 7+1 dataset；只能借鏡統計範式 |
| **R1-m1** (haplotype-specific CNV × GHIR) | 部分 (negative evidence) | A6/A8 對 single-sample ONT；reviewer 想要 LongPhase-S 8 dataset 通用結論，我們只能引「LOH ≠ phasing failure」概念 [I] |
| **R1-m3** (block size distribution) | 部分 (method demo) | A7 我們 HCC1395 PS block N50 943 kb [F]；reviewer 要求 LongPhase-S 自己 7+1 dataset 都報，我們只能示範方法 |
| **R2-#3** (more cancer genome, low rearrangement) | 部分 (negative case study) | 我們 ISM filter 自身 LOSO 也 fail [F]；不是支持 LongPhase-S 而是「同領域 common pitfall」誠實 disclosure |
| **R3-M7** (low VAF P/R + detection 下限) | 部分 (mechanism support) | AF<0.4 是我們 marker filter；但 ISM 不直接報 LongPhase-S 的 detection floor，僅提 caller_af dominance [O] |
| **R3-Sp3** (germline phasing tool head-to-head) | 部分 (mechanism reframe) | 我們未跑 LongPhase mainline (germline-only) vs LongPhase-S (`--somaticMode`) 對照；只有 V5/V6 binary 比較 + `--germline-hp-only` flag toggle |

## §5 不對應條目 — 為何不主動回應

| ID | Skip 理由 |
|---|---|
| **R1-m2** (coverage range 寫進正文) | 純排版；無 ISM leverage |
| **R1-m4** (Fig 6 legend overlap) | 純排版；無 ISM leverage |
| **R2-#1** (cellular vs DNA purity 概念) | LongPhase-S purity model 內部 polynomial 設計問題；ISM HP tag downstream 與此正交 |
| **R2-#2** (train/test overlap 事實錯誤) | LongPhase-S 自行引 Methods 4.3.2 駁回即可；ISM 介入會稀釋論點 |
| **R2-#4** (ASCAT/Battenberg/PURPLE benchmark) | 純 LongPhase-S CNV baseline 對照；ISM 無對應 baseline tool |
| **R3-M2** (PacBio HiFi 跨平台) | 我們手上 7 樣本全 ONT；無 HiFi evidence；硬接會 over-claim |
| **R3-M3** (SNV vs indel mechanism) | ISM 完全不處理 indel；介入無 added value |
| **R3-M4** (filter threshold purity-bin 數值) | LongPhase-S `--somaticMode` filter parameter table；無 ISM 連結 |
| **R3-Sp1** (ClairS vs DeepSomatic 受益不同原因) | 純 caller training data / neural arch 差異 |
| **R3-Sp2** (DeepSomatic / ClairS callset overlap Venn) | 純 caller comparison |
| **R3-Sp4 + R3-Fig** (typos + figure detail) | 純編輯類 |
| **R3 Curiosity-Implicit** (Purity scale 方向 / "Called Somatic Variants" 池) | Figure 排版細節 |

## §6 Provenance

**生成時間**：2026-05-24 (handoff prep for HKU ClairS-TO Luo Lab)
**Git HEAD**：`369e8f9` (branch: `refactor/phase1-safety`)
**Source 文件 last-modified**：
- Reviewer report：2026-05-23（291 行）
- HKU collab HTML v4：2026-05-23 commit `369e8f9` (1167 行)
- Memory files：12d-34d old (verify-against-current-code caveat per system-reminder)

**Evidence tier 標註原則**（per `/scientific-rigor` §2）：
- [F] = **F**irst-hand verified（直接從 source code / BAM / VCF / pysam 量化）
- [O] = **O**bservation（cross-sample 或 single-sample 統計觀察，但機制未完全 elucidated）
- [I] = **I**nference / **I**nterpretation（design proposal / reframe / mechanism 推論）
- [U] = **U**nverified（未實機驗證）

**評論強度分布（v6 corrected）**：⭐⭐⭐⭐⭐ × 2 / ⭐⭐⭐⭐ × 4 / ⭐⭐⭐ × 5 / ⭐⭐ × 5 / ⭐ × 7 = **23 entries total**（含 R3-Curiosity-Implicit）。

**Coverage 分布（v6 corrected）**：✅ 強對應 × 5 / ⚠ 部分對應 × 6 / ❌ skip × **12** = 23 total。

**不過度宣稱原則**：
- 任何 [O]/[I] 統計都以「初步觀察 / 可推論」措辭引用，不寫「證實 / 確認」
- LOH-constrained phasing discovery 已 reframe — Paired ΔNG=+0.787 加註「需獨立 phasing-vs-methylation 驗證」
- HPFineNGroups V6 phase B chr19 conditional POSITIVE 標 caveat「chr19 only = 2.16% genome, BAM=V3F binary」
- Phase 2 Cycle 1 LR filter 已 downgraded to NEGATIVE at sample-level（LOSO -0.00012），不引用作 positive evidence
- ClairS-TO Verdict pilot 內部校準 POSITIVE 但 F1 增益 NEGATIVE — 不主動引用，避免 cross-context confusion

**業界源 footer**：
- Reviewer comment translation methodology: PI-internal rubric (no external citation)
- ISM C++ source: `InterSubMod/src/core/MatrixBuilder.cpp` `DistanceMatrix.cpp` `SignificanceAnalyzer.cpp`
- LongPhase-S source: `/big8_disk/liaoyoyo2001/longphase-s/` v1.0.0 (mainline `twolinin/longphase` v1.7.3 clone)
- ClairS-TO source: GitHub `HKU-BAL/ClairS-TO` `shared/param.py` (`flankingBaseNum=16`)
- SEQC2 truth set: `high-confidence_sSNV_in_HC_regions_v1.2.1.vcf` (chr1+chr8+chr19 6,919 PASS-only TP)
