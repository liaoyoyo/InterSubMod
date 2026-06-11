---
title: ASM × SEQC2 CN ground-truth confound disentanglement (HCC1395 pilot)
date: 2026-06-02
status: in_progress
partial_flag: true   # A pilot — single sample HCC1395 only
task_type: A_pilot
tier: "⭐3 single-sample exploratory (PASS-with-caveats); results filled 2026-06-02"
verdict_summary: "甲基/聚類分群差異主要 mutation-linked；CN/ploidy/error confound 在可測通道上大致排除（H-A/C/D/E CONFIRM, H-B 灰帶）"
anchor_commit: 8a5faa8
sample: HCC1395 paired_full (single sample; only sample with SEQC2 CN/LOH ground-truth)
owner: InterSubMod Research
related:
  - InterSubMod/research/tsg_promoter_asm_reviewer/scripts/18_dual_axis_pivot.py
  - InterSubMod/research/tsg_promoter_asm_reviewer/scripts/24_cluster_eval_core.py
  - InterSubMod/research/tsg_promoter_asm_reviewer/scripts/30_H3_loh_cn_stratification.py
  - memory/project_zar1l_brca2_asm_verification.md
---

# ASM × SEQC2 CN Ground-Truth Confound Disentanglement — HCC1395 Pilot

> **One-line question**: 甲基（ASM Δβ + read-level 聚類）的分群差異，是真的與 somatic mutation 關聯，
> 還是被 copy-number / ploidy / LOH / coverage / error-read / alignment confound 出來的？
> 用 **SEQC2 連續 median CN ground-truth** 反查。

> ⚠ **本檔目前是 PRE-REGISTRATION（跑前鎖定）**。Results / Conclusions 段為空（PENDING workflow）。
> 預測一旦鎖定，達 falsification 條件的 verdict **不可事後改寫**（scientific-rigor §7.1 + CLAUDE.md §1 Hard Gate）。

## 0. Scope & provenance

| 項 | 值 |
|----|----|
| Task type | **A pilot**（subset OK；單樣本 HCC1395）|
| Sample | HCC1395 paired_full（唯一有 SEQC2 CN/LOH ground-truth 的樣本）|
| Outcomes | **O1** = per-site ASM Δβ（HP-axis + ALLELE-axis；dualaxis TSV）+ **O2** = read-level 甲基聚類 blind-ARI（BAM pass）|
| CN ground-truth | SEQC2 Masood 2024（short-read NGS）：`Additional_file_3_cnv_gain_cn_median.txt`（gain CN≥3）+ `_4`（loss CN 0-1）+ `ngs_benchmark_cnvs_gain_loss_loh.bed` |
| BAM | `…/HCC1395/paired_full/…/longphase_s/HCC1395_tagged.bam`（HP + MM/ML 甲基 tag；不複製，indexed read）|
| Anchor commit | `8a5faa8` |

## 1. Pre-registration (Confirmatory) — 預測鎖定

| H | 預測（confirmatory）| 否證條件 (falsifier) | Decision threshold |
|---|------|---------|-----|
| **H-A** | HP-axis ASM Δβ **獨立**於 ground-truth median CN（非劑量效應）| coverage-controlled ρ(\|Δβ\|, median_CN) > 0.2 且同向 p<0.05 | ρ≤0.1 或方向相反 → **非 CN-driven**（支持 ASM 真實）|
| **H-B** | 甲基特徵**無法**預測倍體/CN（超過 coverage baseline）| (methylation→CN class AUC) − (coverage-only AUC) > 0.05 | ΔAUC ≤ 0.02 → **甲基不編碼 CN** |
| **H-C** | HP-axis 與 ALLELE-axis 在 matched loci **一致** | (ALLELE\|Δβ\| − HP\|Δβ\|) 隨 median_CN 上升，interaction p<0.05 | 無 scaling → 兩軸 robust；有 scaling → ALLELE-特有訊號 = allele-dosage artifact |
| **H-D** | ASM 訊號**不被** error/alignment 解釋 | resid(\|Δβ\| \| cov+CN) ~ low-MAPQ/high-NM，ρ > 0.2 p<0.05 | ρ≤0.1 → 技術噪音非主因 |
| **H-E** | read-level 聚類 blind-ARI **不隨** CN/error 系統變化（聚類非結構 artifact）| ARI~median_CN 或 ARI~error 顯著 (p<0.05, \|ρ\|>0.2) | null → 聚類反映真甲基結構而非 CN/對齊 |

**Exploratory（事後分析，標 exploratory，不受上表約束）**：canonical tree 個案（BRCA2 等）細節描述；CN class × TP/FP 交叉；regime 重分解。

## 2. Methods（逐步，對齊 §4 DAG）

DAG：`MUT → METH`（目標因果邊）；confounders `{CN/ploidy, LOH, coverage(n_cpg), error-read, alignment} → METH`；
結構鏈 `CN → alignment難度 → error-read`。**HP-axis 設計上 held-constant CN/ploidy/alignment/region**（同位點同單倍型 germline-read vs somatic-read），故為 de-confounded estimator；ALLELE-axis 不控 baseline ASM。

| Q | 做法 | 主判準 / 對應 H |
|---|------|--------|
| **Q1 方法 re-audit** | O1：抽樣 loci 從 raw Level1（`msa_tmp/*/`）重算 β 比對 dualaxis TSV（驗 MAX-collapse / MIN_N=3 / paired Wilcoxon）；O2：跑「驗尺」PC1 模擬 Δβ=0.3 必 ARI>0.5 & NC1 隨機 split 必 ARI<0.15（script 24 內建）| gate：兩側重現才信下游 |
| **Q2 CN dose-response** | per-site 標 SEQC2 median_CN（連續）；Spearman \|Δβ\|~CN，**分 HP/ALLELE 軸**；coverage-controlled = partial Spearman(given n_cpg) + van Elteren within-cov-stratum；signed Δβ~CN 看方向 | **H-A** |
| **Q3 反向預測倍體** | features {Δβ, \|Δβ\|, baseline β, ARI, sil, n_cpg} → CN class（3-class loss/neutral/gain）AUC；ΔAUC = full − coverage-only；過 `/auc-confound-guard` 三關（within-group OLS / strata / permutation）| **H-B** |
| **Q4 兩軸診斷** | matched loci（同位點有 HP+ALLELE 軸）：(ALLELE\|Δβ\|−HP\|Δβ\|) ~ median_CN interaction（stratified / mixed）| **H-C** |
| **Q5 error/alignment** | BAM pass 取 per-read MAPQ / NM / supplementary → per-site median MAPQ、%MAPQ<20、median NM、%supplementary；residualize \|Δβ\| on cov+CN → 殘差 vs error 特徵 | **H-D** |
| **O2 聚類 × CN/error** | BAM pass per locus blind-ARI（observed-only Hamming, k=2 agglomerative+spectral, PERMANOVA R², silhouette；M8 length-placebo）；ARI ~ median_CN / error | **H-E** |

**O2 聚類 pilot 範圍**（避免硬跑全 39K）：HP-axis sig-ASM（p<0.05，capped ~800 by |Δβ|）+ CN-stratified 抽樣（loss/neutral/gain 各 ~150）+ canonical trees（BRCA2 chr13:32315128, chr20:50267392, chr12:31601630, chr9:42376881, chr16:17774746）。Forest 級 aggregate（Q2/Q3/Q4）用 dualaxis TSV 全量。

## 3. Caveats（誠實邊界）

1. **單樣本 HCC1395** → 所有結論 single-sample exploratory tier；不可宣稱跨樣本一般性。
2. **SEQC2 CN 為 short-read NGS（WGS/WES）建立**，本專案 ONT long-read → 跨平台 caveat；median CN「不應視為 gold-standard」（KB note #4）。
3. **error/alignment（Q5）** 用既有 BAM indexed read（不複製、不重跑 MSA）；`TMPDIR=/big7_disk`。
4. **「倍體」操作化 = 區域 CN**（HCC1395 全基因組 ploidy=2.85 是單一數字，per-site 測的是 regional CN）。

## 4. Data provenance（§8.4）

- O1: `genome_survey_v2/asm_dualaxis_tp.tsv`（51,172 records）+ `asm_dualaxis_fp.tsv`（5,150）
- CN: `/big8_disk/data/HCC1395/SEQC2/CNV/`（v3 命名；Masood 2024, Zenodo v4）
- BAM: longphase_s `HCC1395_tagged.bam`
- 既有 H3 beds: `genome_survey_v2/h3_work/{cnLOH,gainLOH,lossLOH}.bed`
- anchor commit: `8a5faa8`

---

## 5. Results（workflow `wf_8fe183fa-b09` 完成 2026-06-02；8 agents, ~112 min）

> 數據 master table：`master_o1_cn.tsv` = **56,320 records**（= 位點 × 軸，**非位點數**）+ `master_o2_error.tsv`（596 BAM-pass loci）。圖：`figures/cn_confound/q2-q5_*.png`。
>
> ⚠ **Provenance 校正（2026-06-02 用戶 catch，已驗證）**：
> - **唯一位點：TP 30,511 + FP 3,644 = 34,155**；TP records 51,171 = 30,511 位點 × 最多 3 軸（11,739 位點 1 軸 / 16,884 位點 2 軸 / 1,888 位點 3 軸）。
> - **SEQC2 HighConf sSNV truth = 39,447**（`high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz`，全 PASS/HighConf）。**分析的 TP set = 30,490 ClairS-TO TP calls**（`filtered_snv_tp.vcf`，recall 77%；~8,957 truth 為 FN 未進分析）。dualaxis 唯一 TP 位點 30,511 ≈ 30,490。
> - 先前 memory/CURRENT_FOCUS 的「全 39,447 TP」**是把 truth 總數誤當分析 TP set** → 已 ledger correction。FP = ClairS-TO 4,842 FP calls → 3,644 有 ASM test。
> - ⚠ **方法 caveat（pseudoreplication）**：先前 Q2/Q4 統計把同位點的 HP1+HP2 records 當獨立 → p 值偏 anti-conservative；effect 極小（ρ=−0.055）故 verdict 應 robust，但**待 per-position 聚合 re-analysis 確認**（見後續 workflow）。

**資料層（見林）**：HCC1395 是**高度擴增基因組** — CN class gain 43,796 (77.8%) / neutral 11,998 (21.3%) / loss 526 (0.9%)；median_CN 分布 cn3=13,878 cn4=17,484 cn5=6,790…（與 ploidy=2.85 一致）；cnLOH(copy-neutral) 10,342。⚠ gain/neutral 嚴重不平衡、loss 極少 → CN-class 比較需 stratified。

### Q1 方法 re-audit — ✅ GATE PASS
- **O1 β 重算 20/20 完全重現**（max abs diff = **0.0**, tol 0.02）。
- ⚠ provenance：raw Level1 已被 `19_full_tp_asm_batch.sh:63 rm -f` 設計性刪除（disk-safe）；audit **從 MSA binary 確定性重生** 20 個 chr19 HP1 位點（用真 SEQC2 REF/ALT，placeholder allele 會錯分 read → 只 1/20 match；真 allele → 20/20 match 且 n_paired_cpg 完全一致）。scope 窄（~0.04% rows, 單 chrom/axis）。
- **O2 驗尺 valid**：PC1 模擬 Δβ=0.3 → ARI=**0.557**（>0.5 ✓）；NC1 隨機 split → ARI=**−0.005**（<0.15 ✓）。

### Q2 CN dose-response（H-A）— ✅ CONFIRM（非 CN-driven）
- HP-axis **partial Spearman ρ(|Δβ|, median_CN | n_cpg) = −0.055**（p=7.7e-5, n=5142）。**反向 + |ρ|≪0.2**。
- 4 個 coverage strata 全一致負向（ρ∈[−0.073, −0.053]）→ 非 Simpson's paradox。
- integer-CN median |Δβ| ~flat 0.082→0.062（CN1→8，**無單調上升**）。
- 解讀：顯著只因 n 大（power），effect 可忽略；**方向與 dose-response artifact 相反**（artifact 會在高 CN 灌大 |Δβ|=正 ρ）→ 強化「非 CN gain/loss artifact」。

### Q3 反向預測倍體（H-B）— ⚠ INCONCLUSIVE（灰帶 0.02–0.05）
- primary HP-sig（n=5111, neutral vs gain）：**dAUC_max(RF)=0.0463 / logistic=0.0086**。0.02 < 0.046 < 0.05 → 灰帶（falsifier 0.05 未觸發，strict confirm 0.02 RF 未過）。
- **全 absolute AUC < 0.58**（RF full=0.566 / cov-only=0.520 / meth-only=0.547）= 低於 discriminability floor。
- within-cov-stratum RF AUC 0.539–0.571 穩定；meth-only permutation AUC=0.549, **perm p=0.001**（弱但非全 null）。
- o2-augmented block（n=436, +blind_ari/silhouette）dAUC=0.144 → 會 REFUTE，但**棄用**（非 pre-reg 特徵、selection-biased 小子集、overfit）。
- 解讀：甲基帶**弱-但-真**的 CN trace（perm p=0.001），但**實務無用**（AUC<0.58）；唯一未乾淨收斂的假說。

### Q4 兩軸診斷（H-C）— ✅ CONFIRM（ALLELE excess 非 allele-dosage artifact）
- excess = ALLELE|Δβ| − HP|Δβ|（matched loci n=19,937）。OLS excess~CN **slope = −0.00196**（p=4.9e-8 但**負向**，與 falsifier 要求的正向相反；r²=0.0015 ≈ 平）。
- Mann-Whitney gain-vs-neutral excess **p=0.967**（無 class 差）。
- HP-vs-ALLELE |Δβ| 相關：Spearman 0.59 overall，**0.78 within cnLOH**（CN=2 平衡區一致性更高 = dosage artifact 的反面）。

### Q5 error/alignment（H-D）— ✅ CONFIRM-partial
- residualize |Δβ| on cov+CN（r²=0.129，**87% |Δβ| variance 殘餘且與 error 無關**）→ 殘差 vs error proxy：median_nm ρ=−0.008(p=0.85) / mean_nm_per_kb ρ=0.052(p=0.21) / frac_supp ρ=−0.042(p=0.30)。**全 |ρ|≤0.05 ≪0.2**。
- ⚠ **MAPQ blind spot**：median_mapq 恆=60、frac_mapq_lt20 恆=0（零變異）→ low-MAPQ/alignment-difficulty confound **無法測（非 refute）**。結論限縮為「**NM/supplementary error 非主因**」，非「無 error confound」。

### O2 聚類 × CN/error（H-E）— ✅ CONFIRM（弱殘餘）
- blind-ARI vs median_CN ρ=0.116(p=0.005)、vs mean_nm_per_kb ρ=0.119(p=0.004) — 顯著但**<0.2 falsifier bar**；vs median_nm/frac_supp 不顯著。
- **blind-ARI ≫ placebo**（length-split null）：median 0.267 vs 0.015，paired Wilcoxon **p=1.8e-43**（n=338）→ 聚類反映**真甲基結構**，非 read-length collider（collider_flag=false）。
- canonical tree **BRCA2 chr13:32,315,128 blind_ari=0.790**（placebo 0.132, n0/n1=34/21）= 強真 cluster（positive control 過）。

## 6. Verdict per hypothesis（對照 §1 falsifier）

| H | Verdict | falsifier 是否觸發 | 一句話 |
|---|---------|------|--------|
| **H-A** | ✅ **CONFIRM**（supported）| 否（ρ=−0.055 反向且≪0.2）| HP-axis ASM 非 CN dose-response artifact |
| **H-B** | ⚠ **INCONCLUSIVE**（灰帶）| 否（dAUC=0.046<0.05）但 strict confirm 0.02 未過 | 甲基帶弱-真 CN trace 但實務無用（AUC<0.58）|
| **H-C** | ✅ **CONFIRM**（supported）| 否（slope 反向）| 兩軸一致，ALLELE excess 不隨 CN 增 |
| **H-D** | ✅ **CONFIRM-partial** | 否（可測 proxy ρ≤0.05）| NM/supp error 非主因；MAPQ 通道無法測 |
| **H-E** | ✅ **CONFIRM**（弱殘餘）| 否（proxy ρ<0.2）| 聚類=真結構（≫placebo p=1.8e-43），ARI 對 CN 有弱殘餘 |

## 7. Conclusion & next step

**核心結論（single-sample exploratory ⭐3）**：在 HCC1395 + SEQC2 連續 median CN ground-truth 下，**甲基（ASM Δβ）與 read-level 聚類的「分群差異」主要是 mutation-linked，而非 copy-number / ploidy / LOH / error-read / alignment confound 的產物 —— 在所有可測通道上 confound 被大致排除**。機制鑰匙：**HP-axis 設計上把 CN/ploidy/alignment/region held constant**（同位點同單倍型 germline-read vs somatic-read），故為 de-confounded estimator；4 個 confirm 全建立在此軸。

**3 個誠實限制（封鎖 tier 升級）**：
1. **單樣本天花板** — HCC1395 是唯一有 SEQC2 CN truth 的樣本；不可宣稱跨樣本一般性（⭐4+ 需 ≥5/7 樣本，違反）。
2. **H-B 灰帶 + MAPQ 無法測** — 兩個殘餘 confound 風險未完全關閉。
3. **跨平台 caveat** — SEQC2 CN 是 short-read NGS，非 ONT gold-standard。

**與既有結論一致性**：對齊 `project_zar1l_brca2_asm_verification`（ASM 真實但單樣本 ⭐3，需 COLO829）；不與任何 NEGATIVE 衝突。

**Next step（不自動啟動，等用戶）**：
- (a) **COLO829 ground-truth 複製** → H-A/B/C 跨樣本（COLO829 是否有 SEQC2-style CN truth 需查；melanoma 非 SEQC2 主樣本）；
- (b) **H-B 收尾** — 測 CN→coverage/region vs CN→methylation 路徑，把灰帶推出；
- (c) **MAPQ blind spot** — 若要關閉，需 source 有變異的 per-read MAPQ（本 BAM MAPQ 飽和=60）；
- (d) 回 phasing 主軸，本 pilot 當 ASM 真實性的 methods 支撐。

> **Provenance（§8.4）**：anchor commit `8a5faa8`；workflow run `wf_8fe183fa-b09`；scripts `40_cn_annotate.py`–`46_q5_error_align.py`；evaluator synthesis `genome_survey_v2/cn_confound/p4_synthesis.json`；互動 HTML 見同名 `.standalone.html`。

---

## 8. Follow-up（workflow `wf_97d18c3b-42a`，2026-06-02；用戶 catch count + 兩新問題）

> 觸發：用戶發現 56,320 ≠ SEQC2 39,447，並問 (Q-context) 倍體是否不一定影響甲基分群 + 兩具體假設、(Q-measurement) 甲基測量法是否需調整 + subclone。scripts `48`–`52`，evaluator `p4b_synthesis.json`。

### 8.1 Provenance 校正 + pseudoreplication（M1）— ✅ H-A robust
- **records vs positions 釐清**（見 §5 校正段）：56,320 records → 34,154 位點（TP 30,511 + FP 3,643）；truth 39,447 / 分析 TP ~30,490。
- **pseudoreplication 影響極小**：per-record partial ρ=−0.055（n=5142）→ **per-position ρ=−0.039**（n=5056），inflation **1.017×**（只 86/5056 位點同時 HP1+HP2 顯著）。**H-A「ASM 非 CN-driven」HOLDS 且 effect 縮小=更穩**。

### 8.2 倍體是否影響甲基分群（Q-context）— 分兩層，答案精確
- **甲基分群結構（clustering / |Δβ|）= CN-independent ✅**（H-A partial ρ=−0.039；H-E blind≫placebo p=1.8e-43）→ **倍體不驅動甲基分群**。
- **但 caller TP/FP rate = 強烈 CN-context-dependent**（這是兩回事，不可混）：logistic `is_tp ~ LOH×CN+cov` — LOH↓(OR=0.36, p=5e-20)、CN-loss 崩(OR=0.12, p=5e-29)、CN-gain↑(OR=1.57, p=6e-5)、cov↑(OR=1.004/CpG)。最差 cell **LOH.loss TP-rate 0.36**(n=239)、最佳 **nonLOH.gain 0.951**(n=18,321)。
- C1 sub-q：LOH 內 abnormal-CN **不比** cnLOH 更 FP（反向小：cnLOH 18.8% vs abnormal 17.1%，OR=0.892 p=0.009）。

### 8.3 你的兩個具體假設（C2，596 BAM-pass loci）— 都不成立，其一**反轉**
- **H-ctx1「LOH∩異常CN∩clustering→FP」❌ REFUTED + 機制反轉**：cell n=84(FP=9) TP-rate 89.3% vs 93.4%，OR=0.59 p=0.18（ns trend）。**反轉發現：LOH 內 high-clustering 反而 TP-enriched**（OR=2.385, p=0.019；TP blind_ari 中位 0.258 > FP 0.168, MW p=0.0075）→ **FP 集中在 LOW-clustering 的 LOH-gain**。
- **H-ctx2「非LOH∩高cov∩clustering→TP」⚠ UNDERPOWERED**：nonLOH 全域只有 2 個 FP → 結構上無法證（cell 100% TP 是 FP 缺席的 artifact，非證據）。
- ⚠ **不可據此建「high blind_ari = FP filter」**（與直覺相反，且 FP 僅 43/596 嚴重 underpower）。

### 8.4 甲基測量法 scrutiny（M2）— 非全域偏差，但特定位點 lossy
- **全域不偏**：binarized(P≥0.5) vs 連續 β **r=0.998**（Spearman 0.992）；dead-zone(ML 50-200) 僅 12.5%。
- **特定位點 lossy（2 機制）**：① **MAX-collapse 丟掉 5hmC**（corr(db_max, 5hmC)=0.033 正交；5hmC>5mC 佔 12%）② **binarization 壓平 intermediate methylation**（BRCA2 連續 β=0.162 vs binarized 0.072，15% intermediate reads）。
- **建議 targeted（非全域）升級**：報 **分離連續 5mC + 5hmC β**；高 intermediate 位點用連續 β。5mC-only 為 canonical（BRCA2 db_5mC=0.164 vs db_5hmC=0.010）。

### 8.5 Subclone（S）— 候選有，但**未確認**
- per-read 雙峰 43/247 loci(17.4%)，但 **24/43 = HP-split**（germ vs som），**19 個 subclone-like 但 0 個 blind_ari 佐證**（M2 雙峰與 C2 clustering 抽到的 loci 幾乎不重疊）。
- **BRCA2 是 UNIMODAL**（dBIC=−6.6）→ somatic hypomethylation 均勻 = 乾淨 ASM，**非 subclone artifact**（反而強化它是 clean call）。
- **measurement_hides_subclone = YES**（特定位點），但 **subclone claim 卡在第二樣本**（需 COLO829 + M2/C2 共用 loci 重測）。

### 8.6 三問 verdict 表

| 問題 | Verdict | 一句話 |
|------|---------|--------|
| Q-count | ✅ 已校正 + H-A robust | records≠positions；pseudoreplication 1.017× 不影響 H-A |
| Q-context 倍體驅動甲基分群？| ✅ **否**（clustering CN-independent）| 但 caller TP/FP rate 強 CN-dependent（兩層分開）|
| H-ctx1 LOH+異常CN+clust→FP | ❌ REFUTED + **反轉** | LOH 內 high-clust 反而 TP-enriched，FP 在 low-clust |
| H-ctx2 非LOH+高cov+clust→TP | ⚠ UNDERPOWERED | nonLOH 全域僅 2 FP，無法證 |
| Q-measurement 需調整？| ✅ targeted（非全域）| 5hmC 被丟 + intermediate 被壓平 → 分離連續 5mC/5hmC |
| Subclone 存在？| ⚠ 候選有未確認 | 19 候選 0 佐證；BRCA2 unimodal=clean ASM；需 COLO829 |

> **Provenance §8**：workflow `wf_97d18c3b-42a`（6 agents）；scripts `48`–`52`；evaluator `p4b_synthesis.json`；per-position master `master_perpos.tsv`；圖 `figures/cn_confound/{m1,m2,c1,c2,s}_*.png`；互動 HTML `…_followup.standalone.html`。

---

## 9. modkit trial — 測量法升級可行性 + 工具定位（workflow `wf_5f8faad4-a08`，2026-06-02）

> 觸發：用戶要求驗證 ONT **modkit** 是否適用於 §8.4 的 targeted 測量升級（分離連續 5mC+5hmC per-read），並比較 modkit vs MSA vs ISM。modkit v0.6.3 裝於 `/big7_disk/liaoyoyo2001/modkit/`。scripts `54`–`55`；evaluator `modkit_trial/modkit_trial_synthesis.json`。

### 9.1 安裝 + 對照驗證 gate — ✅ PASS（AGREE，乾淨）
- modkit v0.6.3 跑我們 modBAM（1657 reads, 45332 rows, 13 loci）。
- **BRCA2 per-read 5mC Pearson = 1.0**；call-set **完全一致**（shared 13,805 / modkit_only 0 / pysam_only 0）；max diff 0.00195 = ONT 8-bit ML 量化。
- **`?` mode 驗證通過**：0/45332 inferred rows → modkit **不把 skipped/unknown 當 canonical**（這是唯一 load-bearing 風險，已關閉）。
- 操作要點：`mod_qual` 是 **[0,1] float**（非 0-255）；`--include-bed` 需 **BED6**（modkit 拒 BED4）；m/h 預設分開 row（無 `--combine-mods`）。
- → **modkit 是可信賴、可互換的 per-read 5mC extractor**（單樣本）。

### 9.2 GMM 重測（modkit 連續分離 β）— ⚠ NEGATIVE/NULL：沒改變 subclone 結論
- 13 loci（12 候選 + BRCA2）：**9 bimodal、3 LOST（全小 n 邊界）、0 GAINED、0 個 5hmC 新結構**。
- **更乾淨的測量反而更保守**（MAX-collapse 在邊界位點輕微膨脹 BIC）；3 個 lost 全是 n=28–68、dBIC 9.30/8.25/7.40 剛跌破 ≥10 cutoff。
- **分離 5hmC 揭露 0 個隱藏結構** → **直接否證「MAX-collapse 藏了 5hmC 亞克隆」**（5hmC 均勻低 ~0.06-0.10，與 M2 corr=0.033 一致）。
- BRCA2 兩法皆 **UNIMODAL**（modkit 連續 minor weight 0.023 = outlier tail）= 乾淨 ASM。
- → **pysam MAX-collapse 本來就沒扭曲 subclone 圖像**（M2 binarized-vs-continuous r=0.998 已預示）；subclone 卡的是**跨樣本**不是測量法。

### 9.3 modkit vs MSA vs ISM — 互補不取代

| 維度 | modkit (v0.6.3) | MSA | ISM (本專案) |
|------|------|-----|------|
| 角色 | modBAM → per-read 連續抽取 / pileup | **somatic ASM 軸**（HP + ALLELE，tumor+normal）| **per-read subclone 聚類** |
| 5mC/5hmC | **預設分開**（m/h rows）| 都帶，但下游 MAX-collapse **有 bug** | 5mCG-pattern 為主 |
| somatic-aware | ❌（dmr 只是 region-level diff）| ✅ 核心 | partial（TP/FP + HP context）|
| subclone-aware | ❌ | partial | ✅ 唯一專做 |
| `?`-mode | ✅ spec-compliant | 依 modBAM | 依 modBAM |
| 可引用/reviewer 信任 | ✅✅ ONT 官方 | internal | internal |
| 速度 | 快（Rust 多執行緒）| C++ 多執行緒 | C++ O(reads²) 距離較重 |

- **建議架構**：**modkit extract（前端）→ MSA 定 somatic 軸 → ISM 做 subclone 聚類**。modkit 只取代「原始 MM/ML 測量步」+ **順便退役 MSA 的 5mC+5hmC double-row MAX-collapse bug**（−0.054 砍半 artifact 的根因）。

### 9.4 modkit trial verdict 表

| 問 | Verdict | 一句話 |
|----|---------|--------|
| ① 裝+對照 | ✅ PASS（AGREE）| Pearson 1.0、call 完全一致、`?`-mode 正確 |
| ② vs MSA/ISM | 互補 frontend | modkit 抽取（官方可引用）→ MSA 軸 → ISM 聚類 |
| ③ 能否解決 | 能供抽取，但 GMM **無改變** | 測量法本來就夠好；subclone 卡跨樣本非測量 |

### 9.5 結論 + next（不自動啟動）
- **採用 modkit 當 citable extraction frontend**（退役 MSA collapse bug），MSA 留 somatic 軸、ISM 留 subclone 聚類。
- GMM 重測當 **NEGATIVE/NULL 確認** —— 不可僅憑測量法 reopen/re-grade subclone 候選。
- next：(a) 把 modkit-extract 寫進 pipeline manifest（附 cross-val log 作 provenance）；(b) **COLO829 跨樣本** + 報 dBIC bootstrap/CI（3 個 flip 是 n-driven 非生物）。
- ⚠ **count provenance 釐清**：**19** = genome-wide subclone-like（M2，43 bimodal 中）；**12** = s_subclone q3 array 實際持有的 top 候選；**13** = 12 + BRCA2 實測。三者勿混。

> **Provenance §9**：workflow `wf_5f8faad4-a08`（4 agents）；modkit v0.6.3 `/big7_disk/liaoyoyo2001/modkit/`；scripts `54_modkit_extract_crossval.py` + `55_modkit_gmm_retest.py`；輸出 `modkit_trial/{candidates.bed,modkit_extract.tsv,crossval.json,gmm_retest.json,tool_comparison.json,modkit_trial_synthesis.json}`；圖 `figures/cn_confound/modkit_{crossval,gmm_retest}.png`；互動 HTML `…_modkit_trial.standalone.html`。
