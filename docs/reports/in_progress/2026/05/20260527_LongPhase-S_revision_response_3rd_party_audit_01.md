---
date: 2026-05-27
topic: LongPhase-S revision response 3rd-party audit (multi-agent fresh-context cross-validation)
source_files:
  - /big8_disk/liaoyoyo2001/knowledge/raw/paper-notes/longphase-s-reviewer-response-2026-05-23/report.md (291 lines, 5/23 fact-check)
  - /big8_disk/liaoyoyo2001/knowledge/raw/paper-notes/longphase-s-reviewer-response-2026-05-23/report.html
  - /big8_disk/liaoyoyo2001/knowledge/raw/paper-notes/longphase-s-reviewer-response-2026-05-23/supplement.html
  - /big8_disk/liaoyoyo2001/knowledge/raw/paper-notes/longphase-s-reviewer-response-2026-05-23/action-plan.html
  - InterSubMod/docs/concepts/2026/05/20260524_reviewer_coverage_matrix_01.md (v6, 23 entries × evidence tier × ISM coverage stance)
  - /big8_disk/liaoyoyo2001/longphase-s/src/somatic_haplotag/TumorPurityEstimator.cpp:65 (6-coef polynomial source-verified)
agents:
  - R1 audit (aeb8958ec837687a1) — verdict ⭐2/5; Curiosity FAIL unforced error
  - R2 audit (a301565aaccc4657e) — verdict ⭐2/5; LOSO CV missing最致命
  - R3 audit (a1cd735f95e25165a) — verdict ⭐2.5/5; PacBio HiFi + HP3 + multi-subclone 三條 fatal
framework: "Verdict-Pyramid + Cross-Validation Matrix"
purpose: "對 LongPhase-S Ho et al. (HKU 黃耀廷 lab) 對 reviewer 1/2/3 的 revision response 做第三方 audit；cross-validate 5/23 paper-notes 整理 vs 3 fresh-context agent 結論；給 paper-scope GO/PROBE/NO-GO + action items"
audience: PI + HKU collab + 廖子游
git_head: 369e8f9 (refactor/phase1-safety)
rigor_tier: Tier 4 paper-scope (對齊 /scientific-rigor §2-§7 每 claim 標 L1-L5)
overall_verdict: "⭐⭐ 2/5 — Major revision required；當前 rebuttal 結構不足以 secure acceptance"
status: in_progress
---

# LongPhase-S Revision Response 第三方稽核整合報告

> **Tier 4 paper-scope 評估**。3 個 fresh-context audit agent + 5/23 paper-notes fact-check + 5/24 ISM coverage matrix v6 三方 cross-validation。

## §0 TL;DR — Verdict-Pyramid 頂層結論

**整體 verdict ⭐⭐ (2/5) — Major revision required；當前 rebuttal 結構不足以 secure acceptance**。

- **3 agent fresh-context 一致**：R1 ⭐2 / R2 ⭐2 / R3 ⭐2.5 — 與 5/23 report.md 評估「不是全部重寫等級但需補實驗+補 Supp」**基本一致**，但 fresh-context 抓出 **3 條 5/23 沒談到的 vulnerability**（見 §3）。
- **3 個 fatal gap (跨 agent 一致排序)**：
  1. **R3-M2 PacBio HiFi** — COLO829 HiFi 已公開（SMaHT + Zenodo + DeepSomatic CASTLE），作者「promise to add」不可接受。R3 拒絕機率 ~70% [L1: verified public dataset]
  2. **R2-#2 overfitting** — 作者「預估 parameters 不會差太多」+ 沒實際跑 leave-one-cell-line-out CV (LOSO)，被 reviewer 抓的最 vulnerable 句 [L2: industry best practice gap]
  3. **R3-M5 HP3 量化** — 屬編輯部標準補表，response 空白 = 直接拒絕信號 [L1: standard supplementary obligation]
- **ISM 介入策略 (跨 agent 一致)**：**minimal-to-targeted footprint**。
  - 主要切入點 = R1-Curiosity (TSG×methylation IGV) + R3-M5 (HP3 bimodal forward-reference) + R3-M6 (linear depth schema future-work)；
  - R2 全 4 條 ISM 不應主導（fundamental challenge 屬 LongPhase-S 自身義務）；
  - ISM LOSO catastrophe 可在 R2-#2/#3 作 "common pitfall disclosure" 限 1-2 句 limitation footnote。

---

## §1 5/23 paper-notes 整理 vs 現在 author response 對照

> Source: `/big8_disk/liaoyoyo2001/knowledge/raw/paper-notes/longphase-s-reviewer-response-2026-05-23/report.md` (291 lines, 6 共通定錨 fact-checked)

| Reviewer | 5/23 fact-check 評估 | 現在 author response | Δ Verdict |
|---|---|---|---|
| **R1** (Yilei Fu) | 「全部合理且建設性，沒要求重做實驗，主要是補論述+補 Supp 圖表」 | 7 條中 4 條 punt（M1 future / Curiosity 留白 / m1 文字 reframe / m2/m3 標 owner 待補） | 5/23 過於樂觀；fresh-context 升級 Curiosity 為 FAIL ⭐⭐⭐⭐⭐ unforced error |
| **R2** (最嚴厲) | 「主張 #2 事實錯誤直接糾正；其餘合理需補實驗」 | 改名 "tumor DNA fraction" + 提案 chr-subset 拆分 + 1-chr Battenberg + bundled PURPLE | 5/23 對 #2 標「事實錯誤」太硬；fresh-context 提醒：即使 Methods 4.3.2 typo 是 paper 錯，作者「預估參數不差」+ 沒實際跑 LOSO 仍 vulnerable |
| **R3** (使用者導向) | 「全部偏補 Supp + 補 Discussion，沒要求重做核心模型」 | 11+ 條中**僅 3-4 條有實質回應**（PacBio bundled / WhatsHap comparison / DeepSomatic Venn 提及待補） | 5/23 預期 acceptance；fresh-context 升級 PacBio + HP3 + multi-subclone 為 fatal gap |

**對照差異主因**：5/23 是「可達 acceptance 的 path map」；fresh-context audit 是「目前 response 是否已足夠」——後者標 35-40% completion。

---

## §2 三 Agent Cross-Validation Matrix — 逐 reviewer 條目

> 每格 = (Agent Verdict / Evidence Tier L1-L5 / 最關鍵漏洞)
> **L1** = source-verified 量化；**L2** = cross-sample 觀察；**L3** = 機制推測；**L4** = claim only；**L5** = 純設想

### R1 (7 條, agent verdict avg ⭐2)

| ID | R1 agent | 5/23 對比 | 致命度 |
|---|---|---|---|
| **R1-M1** Subclonal/region GHIR | **WEAK / L4** — rename+future punt 不解決量化請求 | 5/23: ✅需補 Discussion+Supp Fig | 中 |
| **R1-M2** Read-level fraction | PASS / L1 — 方法對需補 8-row table | 5/23: ✅容易答 | 低 |
| **R1-m1** Hap-CNV × GHIR | WEAK / L3 — 一段文字 vs reviewer 要 figure | 5/23: ✅補 Supp Fig | 中 |
| **R1-m2** Coverage range | STRONG_PASS / L1 — done | 5/23: trivial | — |
| **R1-m3** Block N50 | PASS / L1 — 標 owner 待補 | 5/23: ✅補 Supp Fig | 低 |
| **R1-m4** Fig 6 重疊 | STRONG_PASS / L1 — done | 5/23: trivial | — |
| **R1-Curiosity** TSG×methylation | **FAIL / L5** — 完全留白 = unforced error | 5/23: 「非強制 IGV 就行」 | **高 (浪費 ISM hero asset)** |

**R1 最弱條目**：Curiosity。ISM 已有 `MatrixBuilder.cpp` + `DistanceMatrix.cpp` + 7-sample ΔNG=+0.787 [L1 量化]，作者**完全沒用**。R1 主動示好的 bonus 不接 = 最廉價失分。

### R2 (4 條, agent verdict avg ⭐2)

| ID | R2 agent | 5/23 對比 | 致命度 |
|---|---|---|---|
| **R2-#1** Cellular vs DNA fraction | PARTIAL / L3 — reframe 自洽但 WGD 對稱假設 hand-wave 無 numerical evidence；training-label systematic bias 沒修 | 5/23: ⚠ 概念對但用錯地方 | 中-高 |
| **R2-#2** Train/test overlap | **PARTIAL / L2** — 拆 chr-subset 不等於 leave-cell-line-out (LOSO)；R²=0.96 in leave-k-out 因 6 cell line 同時出現 train+test 而 over-optimistic | 5/23: ❌「事實錯誤」 | **致命** |
| **R2-#3** Extreme samples | PARTIAL / L2 — 補低 CNA cancer 樣本仍缺；pediatric/heme 缺席 | 5/23: ✅合理但軟性回應 | 中 |
| **R2-#4** ASCAT/Battenberg/PURPLE | PARTIAL / L2 — single-chr Battenberg underpowered；PURPLE long-read 弱應正面提及 (Song 2025 PMC12677682) | 5/23: ASCAT透明✅/Battenberg⚠/PURPLE✅ | 中-高 |

**R2 最弱條目**：#2。作者「**預估 parameters 不會差太多**」是 vulnerable 句，必抓。**單一最高 ROI action = 跑 leave-1-cell-line-out CV (LOSO)** — 同時 address #2 + #3。

**ISM 鏡像警告 [L1]**：`project_phase2_cycle1_global_fp_filter.md` HCC1395 in-distribution 5-fold +0.02236 → LOSO held-out **−0.00012** (100% effect 蒸發, HCC1954 transfer −0.377)。**這支持 R2 擔憂**而非反駁；ISM 應誠實在 Discussion 加 1-2 句 "common pitfall" 而非自誇。

### R3 (11+ 條 + figure, agent verdict avg ⭐2.5)

| ID | R3 agent | 5/23 對比 | 致命度 |
|---|---|---|---|
| **R3-M1** Rescued/removed class | EMPTY_RESPONSE / L4 | 5/23: ✅需 Supp 分類圖 | **高** |
| **R3-M2** PacBio HiFi | **EMPTY/捆綁 R2 / L1 — 已公開可得不能 skip** | 5/23: ✅補 COLO829 HiFi | **致命 (~70% reject)** |
| **R3-M3** SNV vs indel | EMPTY_RESPONSE / L3 | 5/23: ✅補 Discussion 段 | 中 |
| **R3-M4** Filter threshold | EMPTY_RESPONSE / L1 | 5/23: ✅Supp Table | 中 |
| **R3-M5** HP3 fraction | **EMPTY_RESPONSE / L1 標準補表必需** | 5/23: ✅補 Supp Table | **致命 (~60% reject)** |
| **R3-M6** Multi-subclone | **EMPTY_RESPONSE / L4 必補 IGV** | 5/23: ✅IGV 範例 | **致命 (~50% reject)** |
| **R3-M7** Low VAF P/R | EMPTY_RESPONSE / L1 | 5/23: ✅Supp Fig VAF bin | 中-高 |
| **R3-Sp1** ClairS vs DeepSomatic 原因 | 薄弱已提及 / L3 | 5/23: ✅連結 R3-M1 | 低 |
| **R3-Sp2** Callset Venn | 待補 / L1 | 5/23: ✅Supp Table | 中 |
| **R3-Sp3** WhatsHap LOH 對比 | 待補 / L2 — 業界文獻 gap | 5/23: ✅借 LongPhase 原版 | **高** |
| **R3-Sp4** Typos / Figure | EMPTY / L1 — 不修 = unprofessional | 5/23: trivial | 低-中 |

**R3 最弱條目**：M2 (PacBio) + M5 (HP3 量化) + M6 (multi-subclone)。

---

## §3 Fresh-context 抓出的 3 條 5/23 沒談到的 vulnerability（業界文獻佐證）

> 這是 multi-agent fresh-context 主要 added value — 5/23 原稿是「作者立場可如何 spin」，但 fresh-context 用業界文獻找 dock points。

### Vulnerability #1：「tumor DNA fraction」非業界主流 term [L2 文獻]

- **業界 vocabulary survey**：
  - ASCAT (Van Loo Nat Comm 2020) 用 "**Aberrant Cell Fraction (ACF)**" — 定義即 "amount of tumor DNA in sample"，概念**等同**作者新詞但**字面**已被先佔
  - PURPLE (Hartwig) 用 "**purity**" + "ploidy"
  - DeCiFer / PyClone 用 "**cancer cell fraction (CCF)**" per-mutation
  - ctDNA 領域（Foundation Medicine, Adalsteinsson 2017 PMC8165771）**已用 "tumor fraction (TF)"** ← 對作者立場最有利
- **質疑**：為何另立新詞而非對齊 ACF / TF? Reviewer 任一可問。
- **建議**：rebuttal 加 1 句 "this aligns with 'tumor fraction (TF)' in ctDNA literature (Adalsteinsson 2017) and 'ACF' in ASCAT framework" 對齊既有 vocabulary。

### Vulnerability #2：WGD 對稱假設未經 numerical validation [L2 文獻]

- **作者立場**：「WGD 兩 haplotype 同步加倍，HP1/HP2 ratio 不變，GHIR 反映 DNA fraction 而非 ploidy」
- **業界反例 (Bielski et al. Nat Genet 2018; Quinton et al. 2021)**：WGD 常伴 post-WGD chromosome-specific loss → **asymmetric HP retention**
- **HCC1954 即典型 post-WGD breast cancer**（ploidy ~3.9 + extensive post-duplication LOH），HP1:HP2 在 chr17/chr8 已驗證**非對稱**
- **質疑**：作者 hand-wave "都預測 0.8 是對的"，無 ploidy-controlled simulation；R2 大機率會要求 fix-ratio cross-ploidy experiment

### Vulnerability #3：PURPLE long-read 弱點是 motivation 而非 critique [L1 文獻]

- **5/23 評估**：「PURPLE 已能跑 long-read (PacBio 官方 pipeline)，補進來能直接堵 reviewer」
- **fresh-context 補充**：
  - Song et al. *Adv Sci* 2025 PMC12677682 **明確 documented**「all existing CNA tools (incl. PURPLE) performed poorly on long-read」
  - nf-core/pacsomatic 整合 PURPLE+PacBio HiFi 但仍 dev branch
- **策略升級**：作者不應只「補一份 PURPLE」，而應**主動引用** Song 2025 強調 "long-read CNA tool gap is the motivation for LongPhase-S"——把 reviewer #4 的批評**轉成支持立場**。

---

## §4 整合 Action Items — 排序最高 ROI（cross-agent consensus）

| Tier | Action | 對應 reviewer | 工作量 | 影響 |
|---|---|---|---|---|
| **P0 必做** | **Leave-1-cell-line-out CV (LOSO)** — 6 cell line × 5 ratio = 30 points → 6 fold，每 fold hold 1 cell line 全部 5 ratio | R2-#2 + R2-#3 同時 | 2-3 天計算 | **單一最高 ROI**；不做必被 reject |
| **P0 必做** | **COLO829 HiFi pilot** — SMaHT Data Portal / DeepSomatic CASTLE 公開資料，至少 1 chr | R3-M2 | 5-7 天 | 不做 R3 拒絕 ~70% |
| **P0 必做** | **HP3 fraction × TP rate table** — 8 datasets stratified | R3-M5 | 1-2 天 | 標準補表，不做拒絕 ~60% |
| **P1 重要** | **Multi-subclone IGV figure + Discussion limitation** — HP1-2/HP1-3 unsupported design limit | R3-M6 | 2-3 天 | 不做拒絕 ~50%；ISM HKU report §2.3 linear depth schema 可作 future-work proposal |
| **P1 重要** | **Rescued/removed variant stratified table** — LOH × GC × repeat × VAF bin | R3-M1 + R3-M7 | 3-4 天 | 標準 supp 不可省 |
| **P1 重要** | **Ploidy-controlled simulation** — 固定 DNA fraction=0.5 跨 ploidy 不同 cell line，看 polynomial output | R2-#1 | 1 天計算 | 補 WGD 對稱假設證據 |
| **P2 加分** | **R1-Curiosity TSG IGV demo** — ISM `MatrixBuilder.cpp` + `DistanceMatrix.cpp` 實機跑 TP53/CDKN2A/RB1 | R1-Curiosity | 3-5 天 | ISM hero asset；R1 主動示好不接 = unforced error |
| **P2 加分** | **5-chr Battenberg + Song 2025 引用 reframe** | R2-#4 | 3-5 天 | PURPLE long-read gap 轉支持立場 |
| **P3 排版** | typos / figure DPI / Purity scale 方向 / "Uknown" / 藍點 confusion | R3-Sp4 + R3-Fig | < 1 天 | 不修 = unprofessional |
| **P3 速贏** | "tumor DNA fraction" footnote 對齊 ASCAT ACF + ctDNA TF | R1-M1 + R2-#1 | < 1 小時 | 字面 alignment |

---

## §5 ISM 介入策略（跨 agent consensus + 5/24 coverage matrix v6 一致）

| ISM 動作 | 哪 reviewer | 強度 | 5/24 matrix 標 |
|---|---|---|---|
| **R1-Curiosity TSG×methylation IGV demo** | R1 | ⭐⭐⭐⭐⭐ HERO — 必做 | ✅ 強對應 |
| **R3-M5 HP3 bimodal forward-reference** | R3 | ⭐⭐⭐⭐ — Discussion 1 段，~60% reviewer 接受 | ✅ 強對應 (T10 24-chr) |
| **R3-M6 linear depth schema future-work** | R3 | ⭐⭐⭐ — Discussion future-work paragraph | ✅ 強對應 (F4v2) |
| **R3-M1 rescued/removed class single-region precedent** | R3 | ⭐⭐⭐ — F7 chr2:18,086,020 14-site 引用 with caveat | ✅ 強對應 |
| **R1-M1 region-level GHIR signature precedent** | R1 | ⭐⭐⭐ — Inner same-hap NG=2 cross-sample [O] | ✅ 強對應 |
| **R2-#2/#3 LOSO failure disclosure** | R2 | ⭐⭐ — 限 1-2 句 limitation footnote，**不可 over-claim** | ⚠ 部分（matrix 內部不一致已標） |
| **R2-#1 cellular vs DNA fraction** | R2 | 🚫 ISM 不介入 | ❌ skip |
| **R2-#4 ASCAT/PURPLE benchmark** | R2 | 🚫 ISM 不介入 | ❌ skip |
| **R3-M2 PacBio HiFi** | R3 | 🚫 ISM 全 ONT 無對應 evidence | ❌ skip |
| **R3-M3/M4/Sp1/Sp2** | R3 | 🚫 純 caller mechanism | ❌ skip |

**ISM 紅線**：
- ISM 是 separate work，**不能取代** LongPhase-S 端 native analysis
- T10 bimodal / F7 / linear schema 屬 [O]/[F]/[I] tier，措辭限「初步觀察 / 可推論 / future extension」**不寫「證實」**
- Phase 2 Cycle 1 LR filter 已 LOSO 降為 NEGATIVE，**禁止**作正向引用
- ClairS-TO Verdict pilot 內部校準 POSITIVE 但 F1 增益 NEGATIVE — **不主動引用**

---

## §6 業界文獻佐證表

| Claim | Source | Tier |
|---|---|---|
| "tumor DNA fraction" 業界已用 (ctDNA TF) | [Adalsteinsson PMC8165771](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8165771/); [Am J Pathology 2022](https://ajp.amjpathol.org/article/S0002-9440(22)00211-5/fulltext) | L2 |
| ASCAT 用 ACF (Aberrant Cell Fraction) | [Van Loo Nat Comm 2020](https://www.nature.com/articles/s41467-020-17967-y) | L2 |
| PURPLE 在 long-read documented weak | [Song Adv Sci 2025 PMC12677682](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12677682/) | L1 |
| nf-core/pacsomatic 整合 PURPLE+PacBio | [nf-co.re/pacsomatic/dev](https://nf-co.re/pacsomatic/dev/) | L1 |
| ClairS-TO 對 somatic VAF 連續性 + DNA-fraction-aware | [Zheng bioRxiv 2025.03.10.642523](https://www.biorxiv.org/content/10.1101/2025.03.10.642523v1.full) | L2 |
| WGD asymmetric post-duplication losses | [Bielski Nat Genet 2018](https://www.nature.com/articles/s41588-018-0165-1); [Dewhurst PMC6072608](https://pmc.ncbi.nlm.nih.gov/articles/PMC6072608/) | L2 |
| Battenberg = ASCAT 後繼 (subclonal CNA) | [Nik-Zainal Cell 2012](https://www.cell.com/cell/fulltext/S0092-8674(12)00527-2); [Van Loo lab Software](https://www.crick.ac.uk/research/labs/peter-van-loo/software) | L2 |
| Subclonal reconstruction benchmark | [Nat Biotechnology 2024](https://www.nature.com/articles/s41587-024-02250-y) | L1 |
| DeepSomatic CASTLE 6 cell-line tumor-normal benchmark | [Park biorxiv 2024.08.16.608331](https://www.biorxiv.org/content/10.1101/2024.08.16.608331v1); [Nature Biotech 2025 PMC11370364](https://pmc.ncbi.nlm.nih.gov/articles/PMC11370364/) | L1 |
| LRSomatic cross-platform HiFi+ONT | [biorxiv 2026.02.26.707772](https://www.biorxiv.org/content/10.64898/2026.02.26.707772v1.full) | L1 |
| HiFi-somatic-WDL COLO829 Zenodo public | [PacificBiosciences/HiFi-somatic-WDL](https://github.com/PacificBiosciences/HiFi-somatic-WDL) | L1 |
| SMaHT COLO829 FiberSeq HiFi public | [SMaHT Data Portal](https://data.smaht.org/docs/additional-resources/pipeline-docs/long-read_pacbio_hifi) | L1 |
| LOSO / leave-profile-out best practice | [Király J Chemom 2025]; [Ambroise BMC Bioinf 2002](https://link.springer.com/article/10.1186/1471-2105-7-91) | L2 |
| LongPhase 原版 (germline-only) Bioinformatics 2022 | [Lin et al. 38(7):1816](https://academic.oup.com/bioinformatics/article/38/7/1816/6519151) | L1 |
| WhatsHap `--block-cut-sensitivity` low-het cut | [WhatsHap 2.8 docs](https://whatshap.readthedocs.io/en/latest/guide.html) | L1 |
| Distributional bias in LOOCV | [Tougui arxiv 2406.01652](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12662204/) | L2 |
| ABSOLUTE cellular purity definition | [Carter Nat Biotech 2012 PMC4383288](https://pmc.ncbi.nlm.nih.gov/articles/PMC4383288/) | L2 |
| Comprehensive somatic SV benchmark (Illumina/HiFi/ONT) | [biorxiv 2025.09.18.677206](https://www.biorxiv.org/content/10.1101/2025.09.18.677206v1.full) | L1 |
| COLO829 somatic truth set (35,543 SNV) | [Sci Rep srep24607](https://www.nature.com/articles/srep24607) | L1 |
| LongPhase-S preprint | [bioRxiv 2025.11.20.689492](https://www.biorxiv.org/content/10.1101/2025.11.20.689492.full.pdf) | L1 |

---

## §7 Provenance + 推論過程

**評估方法**：
1. **5/23 原稿** (291 lines) fact-checked 6 個共通定錨 — paper Methods 4.3.2 train/test split / `TumorPurityEstimator.cpp:65` 6-coef 多項式 / in silico mixing `samtools view -s` / SNV+4.5% indel+7.1% / HP3 設計 — 均屬 [L1] direct verification
2. **5/24 ISM coverage matrix v6** — 23 entries × {評論強度 ⭐1-5, 嚴重度, 我們對應狀態 ✅/⚠/❌, evidence tier [F/O/I/U]}
3. **Fresh-context 3 agent** — 每 agent 獨立讀 prompt + WebSearch 業界文獻 + ISM coverage matrix
4. **Cross-validate** — 比對 5/23 評估 vs fresh-context；標出 agreement (R2-#2 + R3-M2 + R1-Curiosity gap) + 升級 (5/23 樂觀 → fresh-context 升 vulnerability)

**Provenance**:
- 5/23 原稿: `/big8_disk/liaoyoyo2001/knowledge/raw/paper-notes/longphase-s-reviewer-response-2026-05-23/report.md` (Yao-Ting Huang lab）
- 5/24 整理: `InterSubMod/docs/concepts/2026/05/20260524_reviewer_coverage_matrix_01.md` v6 (23 entries with T10 5/24 升級)
- 3 agent fresh-context transcript:
  - R1 audit: aeb8958ec837687a1 / `tasks/aeb...`
  - R2 audit: a301565aaccc4657e / `tasks/a30...`
  - R3 audit: a1cd735f95e25165a / `tasks/a1c...`
- LongPhase-S source: `/big8_disk/liaoyoyo2001/longphase-s/src/somatic_haplotag/TumorPurityEstimator.cpp:65`
- ISM evidence memories (last verified ages):
  - `project_loh_constrained_phasing_discovery.md` (35d)
  - `project_loh_subclone_af_methylation_positive.md` (35d)
  - `project_hpfinengroups_subclone_marker.md` (16d)
  - `project_phase2_cycle1_global_fp_filter.md` (9d)

**不過度宣稱原則**：
- 任何 [O]/[I] 統計都以「初步觀察 / 可推論」措辭引用，不寫「證實 / 確認」
- LOH-constrained phasing discovery 已 reframe — Paired ΔNG=+0.787 加註「需獨立 phasing-vs-methylation 驗證」
- HPFineNGroups V6 phase B chr19 conditional POSITIVE 標 caveat「chr19 only = 2.16% genome, BAM=V3F binary」
- Phase 2 Cycle 1 LR filter 已 downgraded to NEGATIVE at sample-level（LOSO -0.00012），不引用作 positive evidence

---

## §8 P0 必做 3 項 detailed execution plan

### P0-1 LOSO Cross-Validation (R2-#2 + #3)

**動機**：作者「預估 parameters 不會差太多」是 R2 必抓 vulnerable 句。leave-k-out (k=mixing ratio) 在同 cell line 內 split 不算真 cross-sample；必須跑 leave-1-cell-line-out。

**步驟**：
1. 重訓練 setup：6 cell line × 5 mixing ratio = 30 data points
2. 6 fold，每 fold hold out 1 cell line 全部 5 ratio
3. Fit 6-coef polynomial on remaining 5 cell line × 5 ratio = 25 points
4. Test on held-out 5 points，computed predicted-vs-true DNA fraction
5. Report: per-fold R² + RMSE + bias；overall LOSO R² vs leave-k-out R²=0.96
6. 若 LOSO R² < 0.8：reframe scope to "validated for cell-line panel similar to training set"

**預期結果（cross-validate ISM LOSO experience）**：R² 可能 drop 至 0.6-0.8；HCC1954 (WGD) fold 最 vulnerable，因 ploidy outlier。

**工作量**：2-3 天（含 6-coef refit + result interpretation）

### P0-2 COLO829 HiFi Pilot (R3-M2)

**動機**：COLO829 HiFi 公開資料已可用，作者「promise to add」會被 R3 視為未盡義務。

**Public data sources** (verified 2026-05-27):
- [SMaHT Data Portal](https://data.smaht.org/docs/additional-resources/pipeline-docs/long-read_pacbio_hifi) — COLO829 FiberSeq BAM
- [HiFi-somatic-WDL Zenodo](https://github.com/PacificBiosciences/HiFi-somatic-WDL) — COLO829.30X.SV_region.bam
- [DeepSomatic CASTLE](https://www.biorxiv.org/content/10.1101/2024.08.16.608331v1) — 6 cell-line tumor-normal HiFi benchmark

**步驟**：
1. Download COLO829 HiFi tumor+normal pair（至少 1 chr 如 chr8）
2. 跑 ClairS HiFi mode → 跑 LongPhase-S → 跑 GHIR purity estimation
3. 對比 ONT COLO829_ONT (~33×) vs HiFi (~30×) 的 purity estimate consistency
4. Report: 1 figure (ONT vs HiFi GHIR distribution side-by-side) + 1 table (purity estimate ± CI)

**Caveat**：若 HiFi 與 ONT estimate 不一致 → 需 Discussion 解釋 platform-specific phasing differences

**工作量**：5-7 天（含 download + ClairS HiFi mode tuning）

### P0-3 HP3 Fraction × TP Rate Table (R3-M5)

**動機**：屬編輯部標準補表，response 空白 = 拒絕信號。

**步驟**：
1. 對每個 dataset（8 個）of LongPhase-S 跑完 SomaticHaplotagProcess
2. Count HP1-1 / HP2-1 / HP3 reads from output BAM (htslib aux field)
3. Cross-reference with SEQC2 truth VCF (或 cell-line specific truth set) 計算 HP3 reads 中的 TP fraction
4. Report: 8-row table {dataset, total reads, HP1-1 reads, HP2-1 reads, HP3 reads, HP3 fraction (%), HP3 TP rate (%)}

**ISM forward-reference** (~60% reviewer 接受度，limit Discussion 1 段)：
> "Downstream analysis (InterSubMod follow-up [O]) suggests HP3 may exhibit bimodal behavior — well-phased regions show HP3 enrichment for somatic candidates (Type A, 20/23 chr in HCC1395, TP rate ~90%), while HLA/segdup/sex chromosomes (Type B, 3/23 chr) show HP3 as phasing-failure fallback (TP rate <5%). This interpretation requires further cross-sample validation."

**工作量**：1-2 天（純 count + cross-reference truth VCF）

---

## §9 主要交付建議

3 個落地選項，請選擇：

1. ✅ **存檔為 audit report (本檔已生成)** — `InterSubMod/docs/reports/in_progress/2026/05/20260527_LongPhase-S_revision_response_3rd_party_audit_01.md`
2. **拆 3 個 reviewer-specific rebuttal letter draft (英文)** — 給 HKU 黃耀廷 lab 直接用；3 個 .md 各對應 R1/R2/R3，每個 ≤800 字
3. **生成 single-page HTML dashboard** — Plotly + sticky TOC + 4-axis grading bar chart + L1-L5 evidence ladder，給 PI 1-on-1 用
4. **直接動手 P0 必做 3 項** — 啟動 LOSO CV / COLO829 HiFi pilot / HP3 fraction table 任一

---

**報告完。** 後續可基於 §1-§7 直接展開英文 rebuttal letter，或挑 §8 P0 任一項落地實作。
