---
date: 2026-05-24
audience: PI + HKU-BAL ClairS team (Luo Lab)
status: PI final draft
framework: SCQA (McKinsey/Minto) + Coverage Matrix + 2-section Q list
rigor_tier: Tier 4 (PI-scope, scientific-rigor §2-§7 對齊)
purpose: 對 HKU 説明 InterSubMod 既有甲基/phasing-aware subclone 功能 + 評估能否提升 ClairS-TO caller + 雙 section 問題清單
git_head: 369e8f9 (refactor/phase1-safety)
related:
  - InterSubMod/docs/reports/in_progress/2026/05/20260523_HKU_collab_tag_proposal_01.standalone.html
  - InterSubMod/docs/concepts/2026/05/20260524_reviewer_coverage_matrix_01.md
  - InterSubMod/docs/concepts/2026/05/20260523_HKU_collab_主軸與細節_01.md
  - InterSubMod/docs/concepts/2026/05/20260523_HKU_collab_整理資訊_01.md
  - InterSubMod/research/hku_collaboration/findings_5_23.md
  - /big8_disk/liaoyoyo2001/knowledge/raw/paper-notes/longphase-s-reviewer-response-2026-05-23/report.md
evidence_tiers:
  F: First-hand verified (源碼 / BAM / VCF / pysam 量化)
  O: Observation (cross-sample 或 single-sample 統計觀察，機制未完全 elucidated)
  I: Inference / Interpretation (design proposal / reframe / 推論)
  U: Unverified (未實機驗證)
---

# ClairS-TO HKU Luo Lab Handoff — InterSubMod 甲基/Phasing-aware Subclone 功能説明 + 提升 caller 評估 + 問題清單

> ✅ **T10 全 24 chr 完成（2026-05-24 18:00 update）— 重大 finding: HP3 是 bimodal**
>
> 經 T10 全 24 chr 擴展 scan，原 3-chr 結論「HP3 是 somatic-evidence bucket」修正為**雙峰 (bimodal)**：
> - **Type A**（20/23 chr "well-phased"）：HP3 frac 0.22-0.38% / TP rate 78.6-95.4%（mean 90.4%）= somatic-evidence bucket，~6,900× random background enrichment
> - **Type B**（3 chr "phasing-failure"）：chr6 HLA (9.41% / 1.37%) + chr16 segdup (8.58% / 3.53%) + chrX (21.60% / 0.00% — SEQC2 chrX truth=0) = fallback bucket
> - **Pooled (all 24 chr)**：HP3 1.76% / TP 14.30%（outliers 主導 — 占全 HP3 reads 85.4%）
>
> 對 HKU 的 R3 M5 response 描述須改為「HP3 在 well-phased 區是 somatic-enriched，在 HLA/segdup/sex chr 區是 phasing-failure fallback」雙重定義，避免單一過度宣稱。T11 「見樹也見林」4-layer 整合圖（aggregate / canonical chr1 / extreme outliers chr6+chr16+chrX / well-explained F7 chr2:18,086,020）已產出 `InterSubMod/research/hku_collaboration/figures/T11_4layer_seetree_and_forest.png`。完整 source: `InterSubMod/research/hku_collaboration/findings_5_24/T10_HP3_24chr_findings.md`。
>
> 方法論修正 per memory `feedback_observation_scope_default_comprehensive.md`：完整驗證任務 default scope = 全基因組 + 見樹也見林 4 層。本次因初版繼承 5/23 deliverable A 3-chr scope 未 surfac 提醒，已補強並寫入 feedback memory 避免未來重蹈。

## §1 TL;DR — 3 句結論

- **既有功能準備好對接**：對齊 Reviewer 1 §1.7 已點名 InterSubMod 的延伸方向，我們的 HP 細分樹（linear depth notation `HP1-1-1 / HP1-1-2 / HP1-2-1`）+ ISM read×CpG 甲基聚類補完規則（R/A/X 9-state matrix）可直接 encode multi-subclone 結構進 BAM aux field，C++ 模組（`InterSubMod/src/core/MatrixBuilder.cpp` + `DistanceMatrix.cpp` + `SignificanceAnalyzer.cpp`）已 ready [F]。
- **5 條 reviewer 強對應**：R1-Curiosity (TSG×methylation) / R3-M6 (multi-subclone) / R1-M1 (subclonal GHIR) / R3-M5 (HP3 fraction) / R3-M1 (rescued variant region) — 我們有 L1-L2 實證可主動引用作 LongPhase-S Supp [O+F]，且最強的 R1-Curiosity 是 reviewer 自行 cue 我們延伸研究方向，是最自然對話切入點。
- **誠實 caveat**：我們 ISM filter 自身 LOSO 也 NEGATIVE（HCC1395 in-distribution 5-fold OOF ΔF1 +0.02236 / HCC1395 LOSO held-out ΔF1 **−0.00012** = 100% sample-level circularity bias / 5-sample LOSO mean ΔF1 **−0.00004** / HCC1954 catastrophic transfer **−0.377**）— 不作 positive evidence 而以「同領域 cross-sample failure」誠實 case study 提出，對應 R2 #3 reviewer 對 cell-line 外推 patient 的疑慮 [F negative]。

## §2 Situation — 為什麼現在做這個 handoff？

5/23 HKU LongPhase-S 論文（bioRxiv 2025.11.20.689492）收到三位 reviewer comments（Yilei Fu / R2 / R3），其中 **R1 §1.7 直接 cite「InterSubMod / MethylSomaticAnalysis 的延伸主題」並建議 TP53 / CDKN2A / RB1 IGV demo 進 Supp**（reviewer 報告 line 79-81）。同時 ClairS-TO ssrs 已 source-verified 為 **flankingBaseNum=16（±16 bp window）+ 18 pileup channel + 7 full-aln channel 全無 methylation channel + 不讀 MM/ML tag**（A1 source verification）。5/24 是我們協助 LongPhase-S 強化 rebuttal + 對齊 future ClairS-TO retraining 的窗口；本報告整合 5/22 plan + 5/23 D2 HTML v5 + 5/24 D3 .md 三階段交付。

## §3 Complication — InterSubMod 既有甲基/phasing-aware subclone 功能

我們已具備三項對接 ClairS-TO 的核心能力 — HP tag schema 擴展、ISM 模組 C++ 實現、跨樣本 subclone marker evidence。

### §3.1 HP tag 5 階層 + linear depth notation schema (design proposal [I])

LongPhase-S mainline v1.7.3 `--somaticMode` 已內建 `HP1 / HP2 / HP1-1 / HP2-1 / HP3` 5 階層 string tag（enum `src/haplotag/HaplotagType.h:97-108` HP1_1=5 / H2_1=7）[F]。我們提議的 linear depth schema 擴展：

| 場景 | 觀察組合 | 編碼 | 用途 |
|---|---|---|---|
| A 純 linear | `{RR, AR, AA}` | HP1 → HP1-1 → HP1-1-1 | 兩位點都 +1 somatic linear 累積 |
| B branching 1st-gen | `{RR, AR, RA}` | HP1 → {HP1-1, HP1-2} | 同 1st-gen 起源、兩 mutation 分支 |
| C branching + linear | `{RR, AR, RA, AA}` | HP1 → {HP1-1, HP1-2} → {HP1-1-1, HP1-2-1} | AA 歸屬 ambiguous → ISM 甲基聚類補完 |
| D HP1+HP2 並行 | 兩 lineage 各自 | germline het root → 兩 haplotype 各自帶 A/C 結構 | 全圖 mirror |

**R/A/X 9-state matrix**：兩 SNV pair 觀測（RR / RA / AR / AA + RX / AX / XR / XA + XX = 9，X = read 不 cover 該位點）。**4 缺失補完規則**：`RX→{RR,RA}` / `AX→{AA,AR}` / `XR→{RR,AR}` / `XA→{AA,RA}`，用 ISM read×CpG 甲基相似度聚類 + 鄰近 SNV vote 決定（[I] design proposal，未跑 multi-subclone benchmark）。

### §3.2 ISM read × CpG matrix C++ 模組已實現 [F]

| 模組 | 路徑 | 功能 |
|---|---|---|
| MatrixBuilder | `InterSubMod/src/core/MatrixBuilder.cpp:17-39` | read × CpG matrix builder, MM/ML tag 解析 |
| DistanceMatrix | `InterSubMod/src/core/DistanceMatrix.cpp:310-345` | 6 metric source-verified: **NHD (Normalized Hamming, binary) / L1 (Manhattan, raw) / L2 (Euclidean, raw) / CORR (Pearson, raw) / JACCARD (binary) / BERNOULLI (raw)** — **無 Cosine** |
| SignificanceAnalyzer | `InterSubMod/src/core/SignificanceAnalyzer.cpp` | Fisher 顯著性 |

可實機跑 single TSG region demo（TP53 / CDKN2A / RB1）— 對應 R1-Curiosity reviewer 點名的 follow-up Supp。

### §3.3 Subclone marker evidence 5 source tier

| Source | 觀測 / 數據 | Tier | 路徑 / 參照 |
|---|---|---|---|
| HPFineNGroups chr19 Phase B | 768 regions, NG≥3 TP rate 94.7%, NG_off=5→NG_on=2 max cell TP rate 99.2% | [O] pipeline-dependent | `project_hpfinengroups_subclone_marker.md` (12d old) |
| LOH × NGroups Paired | ΔNG=+0.787 跨 7 樣本 paired mode (p<10⁻⁶⁵), median \|r\|=0.755 | [O] reframed | `project_loh_subclone_af_methylation_positive.md` (31d old) |
| LOH-constrained phasing TO | 6 樣本 NG=2 Inner 93-99% same-hap, TP gap +0.37 | [O] cross-sample | `project_loh_constrained_phasing_discovery.md` (31d old) |
| HKU A2_1 PS set TP 分佈 | 85.9% PS ≥2 TP, n=759, median 5 TP, max 90 | [F] HCC1395 single-sample | `findings_5_23.md` §2 |
| HKU A6 / A7 / A8 | A6 6-zone H=234 p=1.4e-48 但方向反；A7 PS N50 943 kb max 18.3 Mb chr1；A8 24-chr ideogram | [F] HCC1395 single-sample | `findings_5_23.md` §8/§11/§10 |

## §4 Question — 能否提升 ClairS-TO caller？

我們可透過 BAM HP tag → ClairS pileup tensor channel 6 (phasing_info) 把 sub-clone 結構 encode 進 future retraining，但有 5 項預期不能解決的限制需誠實提出。

### §4.1 Channel 對接路徑 [F+I]

ClairS-TO ssrs full-aln tensor channel 6 = `phasing_info`（HP tag）是 LongPhase-S → ClairS 唯一接觸點（A1 source verification [F]）。我們的 HP 細分樹 + ISM 補完 encode 進 BAM aux field（如 `HP:Z:1-1-1` + 補 `MC:Z:<methyl_cluster_id>`）後，future ClairS-TO retraining 可直接 read [I]。

### §4.2 預期 caller 受益 — 5-claim 鎖鏈

| # | 鎖鏈步驟 | 證據 | Tier |
|---|---|---|---|
| 1 | ClairS-TO 視窗 ±16 bp → 僅 2.4% SNV pair 落在 window 內 | A2_2 (`findings_5_23.md` §3): 3,963 pairs 中 96 (2.4%) in-window | [F] |
| 2 | 剩 97.6% pair 的 phasing 訊號透過 BAM HP tag → channel 6 進入 model | A1 ClairS source 驗證 channel 6 = phasing_info | [F+I] |
| 3 | median pair distance 17.4 kb >> ONT read length → 跨 PS chaining 是主來源 | A2_2 median distance 17,421 bp; A7 PS N50 943 kb | [F] |
| 4 | 85.9% PS set 含 ≥2 TP → ISM read-level 甲基聚類有足夠樣本 | A2_1 759 PS sets, 652 (85.9%) ≥2 TP | [F] |
| 5 | HP 細分樹 encode sub-clone 結構 → future ClairS-TO retraining 直接吸收 | linear depth schema design + R/A/X matrix | [I] design |

### §4.3 預期不能解決 — 誠實面

| Reviewer 條目 | 為何 ISM 介入不適當 |
|---|---|
| R2 #1 cellular vs DNA purity | LongPhase-S purity polynomial 內部設計問題，與 ISM HP tag downstream 正交 |
| R2 #4 ASCAT / Battenberg / PURPLE benchmark | 純 LongPhase-S CNV baseline；ISM 無對應 baseline tool |
| R3 M2 PacBio HiFi | 我們手上 7 樣本全 ONT，無 HiFi evidence |
| R3 M3 SNV vs indel | ISM 完全不處理 indel |
| **我們 LR filter cross-sample 失敗** | **HCC1395 in-distribution 5-fold OOF +0.02236 / HCC1395 LOSO held-out −0.00012 (100% effect = circularity bias) / 5-sample LOSO mean −0.00004 (Wilcoxon p=0.125 DIRECTION_NEGATIVE) / HCC1954 −0.377 LOSO catastrophe** — 同 LongPhase-S R2 #3 同類教訓，誠實提作 case study 而非 positive evidence [F negative]。Memory SoT: `project_phase2_cycle1_global_fp_filter.md` DOWNGRADED ⭐⭐ L4 |

## §5 Answer — Reviewer Coverage 對照表（5 強 + 6 部分 + 12 skip = 23 entries）

精簡版引用 `InterSubMod/docs/concepts/2026/05/20260524_reviewer_coverage_matrix_01.md` 完整 23-entry 對照（v6 evaluator audit 後 entry count 修正：5+6+12=23）；5 強對應依強度排序：

| Reviewer 條目 | 強度 | 我們 evidence | 行動建議 |
|---|---|---|---|
| **R1-Curiosity** (TSG×methyl) | ⭐⭐⭐⭐⭐ | Paired ΔNG=+0.787 7/7 + V6 chr19 99.2% + ISM C++ ready | 提供 TP53/CDKN2A/RB1 IGV demo [O+F] |
| **R3-M6** (multi-subclone) | ⭐⭐⭐⭐⭐ | linear depth schema + F4v2 4 scenarios + R/A/X matrix + F7 chr2:18,086,020 14-site | 直接給 HP1-1-1/HP1-1-2 命名 + IGV 範例 [I] design |
| **R1-M1** (subclonal GHIR) | ⭐⭐⭐⭐ | LOH-constrained phasing TO 6-sample + A8 24-chr ideogram | 引 Inner same-hap NG=2 signature + chr8 ~96% chr-level LOH (A8 ideogram) / 99.1% bin-level LOH (A2_3 1kb bins) × HP2-1 dominance [O] |
| **R3-M5** (HP3 fraction) | ⭐⭐⭐⭐⭐ | **T10 (5/24 全 24 chr scan) bimodal HP3**: Type A (20/23 chr) HP3 frac 0.22-0.38% / TP rate 78.6-95.4% mean 90.4% **+ ~6,900× enrichment**；Type B 3 chr (chr6 HLA / chr16 segdup / chrX) HP3 frac 8.6-21.6% / TP rate 0-3.5% = phasing-failure fallback。Pooled HP3 1.76% / TP 14.3%（outliers 主導 85.4%）+ A8 per-chr HP normalized stack | 全 24 chr source-verified [F]：HP3 雙重定義 — well-phased 區 somatic-enriched / HLA-segdup-sex 區 fallback；對 HKU 給 T11 4-layer integrated figure 為 prototype HP3 panel template |
| **R3-M1** (rescued variant region) | ⭐⭐⭐⭐ | F7 chr2:18,086,020 4-Layer Evidence + Inner NG=2 TP rate 0.93 vs Outer 0.55 | 引「LOH × HP2-1 集中 region」為 rescued variant class 範例 [O] |

**部分對應 6 條**：R1-M2 (included-read fraction) / R1-m1 (haplotype-specific CNV × GHIR) / R1-m3 (block size distribution) / R2-#3 (cross-sample LOSO negative case study) / R3-M7 (low VAF mechanism) / R3-Sp3 (germline phasing tool head-to-head)。

**Skip 12 條**：R1-m2 / R1-m4 / R2-#1 / R2-#2 / R2-#4 / R3-M2 / R3-M3 / R3-M4 / R3-Sp1 / R3-Sp2 / R3-Sp4+Fig / R3-Curiosity-Implicit — 純排版 / 純 caller mechanism / LongPhase-S in silico mixing 設計，ISM 介入無 added value。

## §6 雙 Section 問題清單

### §6.1 對 HKU 外部 collaboration ask（14 questions = 5+4+3+3）

**Q-Layer 1：linear depth schema 命名 + retraining 整合（5）**

| ID | 問題 |
|---|---|
| 1.1 | `HP1-1-1 / HP1-1-2 / HP1-2-1` linear depth notation 是否採納？建議 schema 命名 final？ |
| 1.2 | 兩位點外推到 N 位點時的 X 比例上限（IQR 下界 = 2, N=3 可能稀疏）— N 上限建議？ |
| 1.3 | ClairS-TO retraining 加 sub-clone tag — 增獨立 channel（→ 8 channel）vs 共用 channel 6 (phasing_info) 哪個更好？ |
| 1.4 | BAM aux field encoding format（`HP:Z:1-1-1` + 補 `MC:Z:<methyl_cluster_id>`）— round-trip 進 ClairS pileup tensor 是否可行？ |
| 1.5 | HKU 端是否有 sub-clone ground truth dataset 可 benchmark？ |

**Q-Layer 2：ASM evidence 進 ClairS retraining（4）**

| ID | 問題 |
|---|---|
| 2.1 | 是否考慮加 ASM evidence（read-level allele-specific methylation）進 model retraining？ |
| 2.2 | 若加，作為獨立 channel 還是融合進 phasing_info？ |
| 2.3 | ClairS-TO ssrs 既用 HCC1937 / HCC1954 / H1437 / H2009 augmentation，是否考慮把 cell-line ASM 訊號 pretrain 進去？ |
| 2.4 | HP1 vs HP2 read count imbalance ratio 是否在 ClairS-TO 已有內部 metric？（取代 HP-missing 作 LOH 間接指標 — 對應 A6/A2_3 cross-validate） |

**Q-Layer 3：LongPhase-S Rebuttal 對齊（3）**

| ID | 問題 |
|---|---|
| 3.1 | R1-Curiosity TSG demo — 是否願意採用 InterSubMod 跑 TP53 / CDKN2A / RB1 single-region IGV 圖進 Supp？ |
| 3.2 | R3-M6 multi-subclone IGV — 我們 F7 chr2:18,086,020 14-site stratification 是否可作 reference？ |
| 3.3 | R3-M5 HP3 fraction — 我們 A8 24-panel ideogram 可作 prototype HP3 panel — 是否需要我們直接幫做 LongPhase-S 8 dataset？ |

**Q-Layer 4：跨 lab data sharing + future work（3）**

| ID | 問題 |
|---|---|
| 4.1 | LongPhase-S 7+1 dataset BAM 是否可分享給我們做 cross-platform validation？ |
| 4.2 | PacBio HiFi 跨平台 — HKU 是否有 HCC1395 HiFi 或 COLO829 HiFi 可一起做？ |
| 4.3 | 後續 follow-up paper 共同作者 / 致謝 安排？ |

### §6.2 對內 audit list（9 questions = 5+4）

**內 A-Layer 1：ISM 既有 evidence 待補強（5）**

| ID | 問題 |
|---|---|
| A.1 | HP3 fraction per chromosome 完整 scan（T1 script 跑完）— 確認 chr16 / chrX HP3 12% outlier 在 6/7 樣本是否一致？ |
| A.2 | Somatic block length per chr ranking（T2 script 跑完）— 與 A7 N50 943 kb cross-validate？ |
| A.3 | linear depth schema 對 PS set TP 數密度的影響 — N=3 細分後 IQR 下界 2 會降到多少？需重跑 A2_1 |
| A.4 | Inner same-hap NG=2 跨 7 樣本驗證 — 目前 6 樣本 TO + 4 樣本 V6，COLO829 BAM 仍 pending |
| A.5 | HP3 中 TP rate — 我們目前只有 HP3 baseline fraction，未做 HP3-stratified TP rate |

**內 A-Layer 2：production readiness（4）**

| ID | 問題 |
|---|---|
| A.6 | cycle 5+ Path A (phase_block_3d / thread_d pivot) vs Path B (chr8-specific zone gate) vs Path C (low-F1 panel) — 哪個優先？ |
| A.7 | caller_af direction-inconsistent 是 LOSO 災難主因 — 是否有 normalization / per-sample calibration 可解？ |
| A.8 | Phase 2 methyl_filter cycle 1-3 caller_af rank 1 主導，methylation 5th-rank — 是否應該 deprecated 或重新 frame 為「complementary signal」？ |
| A.9 | ISM C++ MatrixBuilder + DistanceMatrix 是否需要 perf optimization 才能 production 跑跨 7 樣本？ |

## §7 Provenance + Caveats

**生成時間**：2026-05-24（5/24 D3 階段；5/22 plan + 5/23 D2 HTML v5 + 5/24 D3 .md）
**Git HEAD**：`369e8f9`（branch `refactor/phase1-safety`）
**Source verification 原則**：所有數字經 source verification（reviewer report line-cited + HTML v5 + coverage matrix + memory cross-ref）。Tier 標記 [F]/[O]/[I]/[U] per `/scientific-rigor` §2。

**不過度宣稱原則**：
- [O]/[I] 統計以「初步觀察 / 可推論」措辭引用，不寫「證實 / 確認」
- LOH-constrained phasing discovery 已 reframe — Paired ΔNG=+0.787 加註「需獨立 phasing-vs-methylation 驗證」
- HPFineNGroups V6 chr19 conditional POSITIVE 標 caveat「chr19 only = 2.16% genome, BAM=V3F binary 非 V5」
- Phase 2 Cycle 1 LR filter 已 downgraded to NEGATIVE at sample-level（HCC1395 LOSO held-out −0.00012 / 5-sample LOSO mean −0.00004 / HCC1954 −0.377），不引用作 positive evidence

**Memory file age caveat**：引用的 5 個 memory 12-34 天舊（`project_loh_subclone_af_methylation_positive.md` 31d / `project_loh_constrained_phasing_discovery.md` 31d / `project_hpfinengroups_subclone_marker.md` 12d / `project_phase2_cycle1_global_fp_filter.md` 3d / `project_clairs_to_verdict_pilot.md` 34d）— 引用前已 verify against current code per `.claude/CLAUDE.md` §9 reference manual rule。

**Coverage 分布（v6 audit-corrected）**：✅ 強對應 ×5 / ⚠ 部分對應 ×6 / ❌ skip ×12 = **23 entries total**（17 substantive ≥⭐2；含 R3-Curiosity-Implicit 編輯類）。

**業界源 footer**：
- SCQA framework: Barbara Minto, *The Pyramid Principle*（McKinsey）
- Coverage Matrix methodology: PI-internal rubric
- ISM C++ source: `InterSubMod/src/core/MatrixBuilder.cpp` `DistanceMatrix.cpp` `SignificanceAnalyzer.cpp`
- LongPhase-S source: `/big8_disk/liaoyoyo2001/longphase-s/` v1.0.0 clone of mainline `twolinin/longphase` v1.7.3
- ClairS-TO source: GitHub `HKU-BAL/ClairS-TO` `shared/param.py` (`flankingBaseNum=16`)
- SEQC2 truth set: `high-confidence_sSNV_in_HC_regions_v1.2.1.vcf` (chr1+chr8+chr19 6,919 PASS-only TP)
