# Prior Art Notes — V3F/V5/V6 ISM × LOH × HP × CN 後分析

<!--
建立時間: 2026-05-15
目標: prior-art 差異化分析（v0.3 計畫 §Agent D 任務）
範圍: TumorLens / ROCIT / SGZ / Wakhan / SAVANA 全文閱讀 + 與本研究差異化
KB 查詢結果: 上述 5 個工具均無對應 doc（/big8_disk/liaoyoyo2001/Knowledge/05_tools/ 內僅有 longphase / longphase-to / longphase-s / variant-callers / methyl-somatic-analysis / intersubmod / tools-overview）
取材限制: medRxiv/bioRxiv full-text 全部 403；TumorLens/Wakhan 用 abstract + GitHub README + WebSearch 摘要；ROCIT 用 PMC mirror；SGZ 用 PLoS 全文；SAVANA 用 PMC mirror
-->

## §1 Prior Art 摘要表

| 工具 | 軸（features） | Read vs Variant | TP/FP 評估方法 | Accuracy claim | 主要限制 |
|------|---------------|----------------|---------------|---------------|---------|
| **TumorLens** (medRxiv 2026-03-19) | SNV + indel + SV + CNV + LoH + CpG methylation + HLA；以 sample-level **purity-aware** modeling（adapted Spectre）joint output；每 haplotype split 做 allele-specific methylation 4×4 panel | **Sample-level / region-level**（10 kb methylation windows），非 per-read TP/FP | 與 truth set 對 SV recall / LoH count / variant retention 比 | 90% 純度: SV recall 86%（+10% disjoint）、17 LoH regions；70%: 15 LoH；50% 純度: 82.7% variant retention | Sample-level 評估，未發表 per-read TP/FP；HLA 為核心場景；CNV 模型 = Spectre adaptation 非 first-principles |
| **ROCIT** (bioRxiv 2026-03-05, PMC12991090) | Per-read CpG methylation prob + read-internal CpG 位置 + sample methylation distribution + 84 cell-type reference atlas | **Per-read classification**（transformer encoder, 384-dim × 6-head × 3-block）；training labels 從 LOH 失 allele reads + somatic SNV同 hap reads 推得 | Mean test AUC、Pearson 與 CN-expected proportion；perturbation Δ | **Mean AUC 0.933**（test sets）；read length 15 kb avg 0.92；observed vs expected tumor read fraction Pearson 0.947 (p=8.1e-11) | Per-sample 訓練（cross-sample AUC ↓0.22）；只用 methylation feature，**未做 LOH × HP × CN 三軸 TP/FP grid characterization**；需 PacBio HiFi + Jasmine 5mC；training 用 LOH/CN/HP 作 label generator，**非 grid axis** |
| **SGZ** (Sun et al. PLoS CB 2018, doi:10.1371/journal.pcbi.1005965) | tumor purity p × ploidy Ψ × per-segment copy number C_i × minor allele M_i → expected AF；MCMC Gibbs + grid search 找 (p, Ψ, {C_i, M_i}) | **Variant-level**，short-read FoundationOne panel | Binomial hypothesis test (α=0.01) vs expected germline/somatic AF；binary call rate + accuracy | 85% call rate；somatic 95-99% accuracy；germline 97-99%（vs basic AF-only baseline 67-95% / 41-87%）；無 F1/AUC | 需 ≥10% normal admixture；≥20% purity 才能 zygosity call；**無 phasing/haplotype 軸**；無 methylation；無 long-read；無 per-read TP/FP characterization |
| **Wakhan** (medRxiv 2025-12-15) | phased VCF + tumor BAM → haplotype-specific CN segments；用 HP CN 差異反向修 phase-switch；purity-ploidy grid search | **Segment-level**（haplotype-specific CN）；無 per-variant 或 per-read TP/FP | 與 Illumina 短讀 CN 一致性 + smaller CNA 精確度（無正式 F1/AUC） | "consistent" with Illumina；小 CNA precision 優於短讀；purity range default 0.5-1.0 / ploidy 1.5-5.5 | 不做 SNV/methylation/per-read；輸入需已 phase 的 het VCF；CN segment-level 為終產物，**無 TP/FP TP cell 概念** |
| **SAVANA** (Nat Methods 2025, doi:10.1038/s41592-025-02708-0, PMC12240814) | somatic SV/SCNA from long reads；可選 phased input 達 single-haplotype 解析；用 BAF 在 LOH 區估 purity；**無 methylation** | **Read-level**（scan SV-supporting alignments + read-backed phasing 驗證），但聚合到 SV/CNA **variant call**；不在 SNV 軸 | ROC AUC（vs 99 tumor-normal pairs）；vs 二/三名工具 specificity 比較 | SV mean AUC **0.98** (0.97-0.98)；specificity **13× 與 82×** 高於次佳；86% recall vs Illumina | 不做 SNV（與本研究主軸正交）；不做 methylation；對 intratumor heterogeneity SV 偵測敏感 |

## §2 與本研究的差異化矩陣

| 維度 | TumorLens | ROCIT | SGZ | Wakhan | SAVANA | **本研究 (V3F/V5/V6 後分析)** |
|------|-----------|-------|-----|--------|--------|----------------------------|
| 分析層級 | sample / 10kb region | **per-read** | variant | segment | variant (SV/CNA) | **per-read & per-region** TSV → 50-cell grid |
| Variant type | SNV+indel+SV+CNV+LoH | (none, 純 methylation classifier) | SNV | CNV | SV+CNA | **SNV** (ClairS-TO somatic) |
| LOH 軸用法 | LoH 計數 (detect 為 output) | **training label generator** | 隱含於 CN minor allele M | output segment | purity 估計 input | **post-hoc TP/FP grid 軸**（Inner/Outer × HP bucket × CN bin）|
| HP/phasing 軸 | per-hap methylation split (4×4 panel) | training 用 same-hap-as-mutation | **無** | 用於 CN segment 估計 | optional input | **HP 4-bucket × 3 longphase variant** (V3F/V5/V6 比較) |
| CN/coverage 軸 | Spectre purity-aware adaptation | 用作 expected tumor read frac | 4-axis grid (p, Ψ, C, M) | output | output | **Coverage_Multiple 5-bin KDE-corrected** (cov_loss / normal / elevated / gain / high_gain) |
| Methylation | 10kb window allele-specific | per-CpG per-read prob | **無** | **無** | **無** | post-hoc TSV (master.tsv ISM output)，輔助 covariate（不在主 grid） |
| 三個 phasing variant 比較 | – | – | – | – | – | **V3F vs V5 vs V6**（priority-bug fix 兩段貢獻分解）|
| TP/FP 評估口徑 | recall / retention vs truth | per-read AUC | binary call accuracy | CN concordance | SV AUC / specificity | **per-cell TP rate + FP rate + log-odds + Wilson 95% CI + Fisher p + within-group OLS confound guard** |
| Power-aware stratification | – | – | – | – | – | **powered (n≥50) / marginal (30-49) / underpowered (<30)** 三層 |
| Confound guard | (none stated) | perturbation cluster test | binomial null only | – | replication | **within-group OLS on NumReads + caller_af + permutation null + L3 AF-bin cross-check + spatial autocorr guard** |

## §3 為何本研究仍有新貢獻（5 點）

1. **三向 phasing variant 比較唯一性**：ROCIT/SAVANA/Wakhan 都假設 phasing 為 input（黑盒），TumorLens 走 ClairS+longphase pipeline；**無任何 prior art 在同一資料上對比 V3F (PON-only) → V5 (Layer 1.5) → V6 (germline-absent revert) 三個 haplotag 版本的下游 TP/FP signature**。本研究是 longphase priority-bug fix 對 read-level marker 影響量化的首例。

2. **Post-hoc grid characterization, 非 caller / classifier 設計**：所有 prior art（TumorLens/ROCIT/SGZ/Wakhan/SAVANA）都是**新工具 / 新模型**，自帶 truth-set 評估。本研究**不發新工具**，是在已有 ClairS-TO 輸出之上做 LOH × HP × CN 三軸 TP/FP 解剖 — 屬於 caller behavior characterization 的後分析學科，與 prior art 屬性正交（complementary not competitive）。

3. **Coverage_Multiple KDE-corrected 與 LOH × HP 交叉首次**：Wakhan/TumorLens 把 CN 當 output，SGZ 把 CN 當 expected-AF input；**沒有任何 prior art 把 KDE-corrected per-read coverage 當 grid 軸與 LOH inner/outer × HP 4-bucket 三軸交叉做 cell-level TP/FP rate**。Coverage_Multiple r=0.831 vs SEQC2 是本專案獨有 binary fix（commits 374fad4 + 12d9b3e）。

4. **Power-stratified cell-level grid + 多層 confound guard**：ROCIT 做 perturbation cluster check；SAVANA 做 replication；TumorLens 做 purity dilution；**沒有 prior art 提供 50-cell × within-group OLS (NumReads + caller_af) × permutation null × L3 AF-bin cross-check × spatial autocorr 五層 confound guard 的方法骨架**（這也是本專案吃過 L2 collider bias / pooled OLS / spatial autocorr 三次教訓後沉澱出來的，memory feedback 三條對應）。

5. **HCC1395 chr8 hotspot 機制泛化測試**：LOH+HPMergedSig 7.4× FP enrichment 87.5% 來自 chr8（已知樣本特異）；本研究 Z-CHR8 vs Z-OCH vs Z-GL vs Z-AUTO 四 zone 對照在 prior art 完全未涵蓋 — TumorLens/SAVANA 都是 cohort-level recall/precision，未做 zone-by-chr 機制歸因。

## §4 風險評估 — 若 prior art 全文涵蓋了我們的 50-cell grid

**最壞情況設想**：TumorLens 或 ROCIT 全文（無法取得 — 403）內藏一個 supplementary figure 做了 LOH × HP × CN 三軸 cell grid。

**本研究差異化基礎仍站得住的 4 個理由**：

1. **三 phasing variant 比較不可替代**：即使他們做了 grid，也是用**一個** haplotag pipeline 的結果（TumorLens 用自家流程；ROCIT 用 HiPhase）。**V3F → V5 → V6 兩段差量分解（Layer 1.5 加入 vs priority-bug fix）只能在我們的 phaseC_genome_three_way 12 個 ISM run 上做**。
2. **ClairS-TO 特定 caller behavior**：所有 prior art 用各自 caller（TumorLens 自家、ROCIT 用 SAGE、SGZ 用 FoundationOne pipeline、Wakhan/SAVANA 不做 SNV）。**ClairS-TO 在 LOH × HP × CN cell 內的 TP/FP signature 是 caller-specific，沒人會剛好做過**。
3. **HCC1395 SEQC2 truth 加 KDE-corrected coverage**：我們的 Coverage_Multiple 是用全 7 樣本共建 KDE baseline + r=0.831 vs SEQC2 校準，prior art 的 CN 用各自 caller，**數值口徑不可比**。
4. **Characterization-only 邊界明示**：我們**不評估 filter/ΔF1**（Plan §5 + memory `feedback_paired_f1_filter_abandoned.md`）；prior art 全都做 caller/classifier benchmarking（含 F1/AUC/recall claim）。**研究問題本身不同** — 我們問「TP/FP 在這 50 個 cell 怎麼分布」，他們問「我的工具準不準」。

**若風險真實發生（prior art 隱含同 grid）**：本研究仍可以「ClairS-TO + longphase three-variant + KDE-corrected coverage 三項組合是首次」立足；最壞需在報告 §Related Work 加註「TumorLens supplementary fig X 做過類似 grid，但用 X caller 與 X coverage」。

## §5 引用清單

| 工具 / 文獻 | DOI / URL | 取得狀態 |
|------------|----------|---------|
| TumorLens (medRxiv 2026-03-19) | doi:10.64898/2026.03.18.26348569; https://www.medrxiv.org/content/10.64898/2026.03.18.26348569v1 (full **403**) | WebSearch 摘要 + abstract |
| ROCIT (bioRxiv 2026-03-05) | doi:10.64898/2026.03.03.709085; https://www.biorxiv.org/content/10.64898/2026.03.03.709085v1 (**403**) → mirror **PMC12991090** ✅ | PMC 全文成功 |
| SGZ (Sun et al. PLoS Comp Biol 2018) | doi:10.1371/journal.pcbi.1005965 ✅ | PLoS 全文成功 |
| Wakhan (KolmogorovLab, medRxiv 2025-12-15) | doi:10.64898/2025.12.11.25342098 (full **403**) | GitHub README + WebSearch 摘要 |
| SAVANA (Nat Methods 2025) | doi:10.1038/s41592-025-02708-0; PMC12240814 ✅ | PMC 全文成功 |
| **KB 內既有相關** | `/big8_disk/liaoyoyo2001/Knowledge/05_tools/{longphase, longphase-to, longphase-s, variant-callers, methyl-somatic-analysis}.md` | 已有，但**無上述 5 tool 對應 doc** |
| 本專案 V6 驗證 | `research/paired_priority_bug_audit/09_V6_caller_F1_verification.md` | local |
| Thread D 主報告 | `InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md` | local |
| Self-Phasing 觀察整合報告（V5 缺陷 §8.6） | `InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md` | local |

## §6 KB 補建建議（提案，不執行）

KB `/big8_disk/liaoyoyo2001/Knowledge/05_tools/` 目前缺以下 5 個 tool doc。建議補建（**待用戶決定**，本任務不寫）：

| 建議檔案 | 內容大綱 | 優先序 |
|---------|---------|--------|
| `05_tools/tumorlens.md` | (1) 工具定位（long-read tumor-only joint SNV/SV/CNV/LoH/methyl/HLA）；(2) Spectre 改編的 purity-aware CNV 演算法；(3) 4×4 allele-specific methylation panel 設計；(4) 50%/70%/90% 純度 benchmark 數字；(5) **明確標註不做 per-read TP/FP**；(6) 與本專案差異化句子（避免未來 agent 從 TumorLens summary 推論我們的 read-level 結論） | **高** — 與本專案最接近 |
| `05_tools/rocit.md` | (1) per-read methylation classifier（transformer 384×6×3）；(2) 4 種 feature 來源；(3) AUC 0.933 / Pearson 0.947 數字；(4) **明確標註只用 methylation 軸；LOH/CN/HP 是 training label generator 非 grid axis**；(5) per-sample training 限制 (cross-sample AUC↓0.22)；(6) 與本專案差異化句（我們做 grid characterization, 不是 per-read classifier） | **高** — 容易混淆「也做 read-level」 |
| `05_tools/sgz.md` | (1) FoundationOne 短讀 tumor-only somatic/germline 二分；(2) 4-axis (p, Ψ, C, M) MCMC + grid；(3) 95-99% accuracy；(4) **無 phasing/methylation/long-read**；(5) 與本專案差異化（我們的 LOH × HP × CN grid 與 SGZ 4-axis 都不同口徑） | 中 — 歷史錨點 |
| `05_tools/wakhan.md` | (1) long-read haplotype-specific CN segment caller；(2) 用 HP CN 差異反向修 phase-switch；(3) GitHub default purity/ploidy range；(4) **不做 SNV/methylation/per-read TP/FP**；(5) 與 V6 BAM 流程關係（Wakhan 可選但 KDE-only 即足夠） | 中 |
| `05_tools/savana.md` | (1) long-read somatic SV/SCNA；(2) AUC 0.98 / specificity 13×/82×；(3) **不做 SNV/methylation**；(4) 與本專案 SNV 軸正交；(5) 提醒 ClairS-TO callset 的 SV 不在 SAVANA 範圍 | 低 — 與 SNV 軸正交 |

建議補建順序：tumorlens.md → rocit.md（最容易與本研究混淆）→ 其他三個。每檔 ≤300 行；可放 `05_tools/{name}.md` 與 `05_tools/en/{name}.md` 雙語對照（沿用既有命名）。

---

**Agent D 任務狀態**：prior art 全文閱讀完成（5/5 工具；2 個透過 PMC mirror、1 個 PLoS 全文、2 個用 abstract/GitHub + WebSearch 摘要因 medRxiv 403）；KB 補建檔案清單已列；本研究差異化基礎站得住（即使最壞情況）。
