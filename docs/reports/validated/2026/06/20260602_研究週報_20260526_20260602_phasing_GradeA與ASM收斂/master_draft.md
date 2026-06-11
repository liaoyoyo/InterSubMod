<!--
建立時間: 2026-06-02
狀態: validated (週報母稿 master_draft)
報告類型: weekly_master_draft
受眾: PI
framework: Multi-Thread-Narrative (weekly-report thin wrapper)
涵蓋期間: 2026-05-26 ~ 2026-06-02
git_HEAD: 37e1511
data_sources:
  - research/tpfp_loh_af_kde_discrimination/data/obs18b_wilcoxon_ng2_gap_n7.json
  - docs/experiments/in_progress/2026/05/20260530_LOH_phasing_n7_cross_sample_Grade_A.md
  - docs/experiments/in_progress/2026/05/20260529_ISM_ZAR1L_BRCA2_ASM_verification_01.md
  - docs/experiments/in_progress/2026/05/20260531_ISM_aux_tag_observation_funnel_01.md
  - docs/experiments/in_progress/2026/06/20260602_ASM_CN_confound_disentanglement_pilot_01.md
  - docs/experiments/in_progress/2026/06/20260602_ASM_cluster_characterization_regimeA_LOH_diagnostic_01.md
  - docs/reports/in_progress/2026/05/20260527_HCC1395_TSG_promoter_ASM_reviewer_response_01.md
  - docs/reports/in_progress/2026/05/20260527_LongPhase-S_revision_response_3rd_party_audit_01.md
  - docs/data_specs/20260601_cleanup_audit/FINAL_CLEANUP_LEDGER.md
  - research/autoresearch/evidence_ledger.jsonl
  - docs/CURRENT_FOCUS.md
provenance_note: 全數字經 2026-06-02 gather→verify workflow (wf_02d8afb0-e13, 12 agents) 雙階段對抗 grep 驗證；116/117 通過，唯一 grep-fail 'FDR 3,834(7.5%)' 與衍生比例 '(0.61%)' 已剔除不引用。
-->
<!-- provenance-verified: 全 metric 數字經 gather→verify workflow wf_02d8afb0-e13 雙階段 grep 對抗驗證（每執行緒 gatherer 萃取 + 獨立 verifier 回源頭檔 grep）；source 路徑見 frontmatter data_sources。落地 CLAUDE.md §13 — 寫此檔的 Write 與產數字的分析不同 batch。 -->

# 研究週報母稿 — 2026-05-26 ~ 2026-06-02

> **主軸（≤30 字）**：本週把論文主軸 phasing 推到 Grade A 門檻，並把 ASM 支撐軸全面自我證偽、收斂成 reviewer-proof 的 supplement-null。
> **主線類型**：`進展:問題` 混合 — phasing 進展型（主）+ ASM 收斂/釘死邊界型（sub）。
> **一句話重點**：這是「**收斂與釘死邊界**」的一週，不是「擴張」的一週。

---

## §0 Highlights / TL;DR

⭐⭐⭐ **This Week's Verdict**
phasing 主軸 n=6→**n=7 達 Grade A 要求 #1**（Wilcoxon W=28, p=0.0078125）；ASM 支撐軸經 6/02 CN-confound pilot + 5/31 六-agent 紅隊全面自我證偽後，定位收斂為 **real but non-directional + non-discriminative + CN-confound-排除**的 characterization + supplement-null。

**Top Findings**
- ⭐⭐⭐ **主軸 phasing 升級**：NG=2 Inner same-HP1 vs Outer cross-het 的 TP-rate gap **7/7 樣本一致正向**，加入 COLO829 後 W=28 **p=0.0078125**（exact = 1/2⁷ 理論最小）。[src: obs18b_wilcoxon_ng2_gap_n7.json]
- ⭐⭐⭐ **ASM confound 大致排除**：6/02 pilot 用 SEQC2 CN ground-truth 反查，H-A/C/D/E **CONFIRM**、H-B 灰帶 → 甲基分群差異主要 mutation-linked，非 CN/ploidy/error/alignment artifact。[src: 20260602_ASM_CN_confound_disentanglement_pilot_01.md]
- ⭐⭐ **ASM 不能當 filter（釘死）**：連續 |Δβ| AUC=**0.5049**（perm-p=0.496 中性）；strong-ASM 反在 FP 富集 **5×**（Fisher OR=0.194, p=1.8e-28）= 回歸極值 artifact 非判別力。[src: 20260531_ISM_aux_tag_observation_funnel_01.md]
- ⭐⭐ **regime-A 乾淨子集有真訊號**：collider battery 後 **15/23** 存活，blind-ARI **0.229 vs germline-het 0.038**（Cliff δ=0.368, p=2.3e-4）。[src: 20260602_ASM_cluster_characterization_regimeA_LOH_diagnostic_01.md]

**Top Asks（需 PI 判斷，[U]）**
1. [U] phasing 是否花 ~25-50hr C++ flag-on 全 7-sample 負控，把 Grade B+ 補成 full Grade A？**（本週決議：先寫作，reviewer 要求再補）**
2. [U] ASM ⭐4 跨樣本（COLO829）維持 off-goal（不為 phasing 論文做）？
3. [U] G6 framework draft 是否本週啟動？

**Decisive Next Step**：phasing 主軸證據已足以撐論文核心 → **先寫 G6 framework draft 把 B+ 證據定稿**，full-Grade-A 負控（長計算）不阻塞寫作，留待 reviewer 要求時補。

---

## §1 終極目標 G6 — 三軸合一論文（進度地圖）

**G6 = read-level epigenetic & haplotype characterization 論文**（用戶 2026-05-29 確認 pivot，characterization 而非 variant filter；產出 76-task DAG / 6 final goals / 11 blockers / 27 guardrails，adversarial-verify PASS）。[src: CURRENT_FOCUS.md / project_goal_landscape memory]

```
G6 論文
├─ 主軸   ⭐3 Grade B+→A    LOH-constrained phasing  ── 本週 n=7 W=28 p=0.0078 ✅ 達 Grade A 要求#1
│                                                        （剩 #2：全 7-sample flag-on 負控 ~25-50hr C++）
├─ 支撐   ⭐3               ZAR1L/BRCA2 ASM           ── 本週自我證偽 → supplement-null
│                                                        （real / non-directional / non-discriminative / CN-confound 排除）
└─ methods ⭐⭐⭐⭐ L2       LOSO sample-level NEGATIVE ── 已是核心 claim（誠實方法學）
   ✗ G5 甲基→F1 filter   ⭐2 L4 DEAD（100% circularity；勿再開）
```

**原始 5 目標 → G6 映射**：目標2（clone 結構）→ phasing 主軸；目標1（per-CpG 甲基評分）→ ASM 支撐軸；目標3（二次打擊順序）暫緩（依賴目標1，P4 pilot region-level 無 two-hit 訊號）；目標4（TO normal 補強）off-goal；目標5（evidence 整合）= G6 本身。[src: project_research_vision_five_goals memory，整理判讀]

---

## §2 三根支柱：假設 → 方法 → 結果（Layer 2 證據卡）

### §2.1 支柱一（主軸）— LOH-constrained phasing｜LIVE ⭐3 Grade B+→A ✅[F]

**One-line Verdict**：**NG=2 Inner same-HP1 的 TP-rate 顯著高於 Outer cross-het，方向跨 7 樣本一致，cross-sample 軸達 Grade-A strength。**

| 項 | 內容 |
|---|---|
| 假設 [F] | NG=2 區 Inner same-HP1 TP-rate > Outer cross-het；方向跨樣本一致（TO mode，純 LongPhase + LOH.bed 可重現，不需甲基）|
| 方法 | 每樣本 gap = Inner − Outer TP rate（Extreme-AF 子集）；Wilcoxon signed-rank exact one-sided greater + bootstrap CI（seed 20260423）；n=6 → **加 COLO829 達 n=7** |
| 結果 | **7/7 全正向**；W=**28**；p=**0.0078125**；median gap **0.345**；bootstrap 95% CI **[0.104, 0.459]**（從 n=6 W=21 p=0.015625 6/6 強化）|

**見林 — per-sample gap（canonical 表）**
| 樣本 | gap | 備註 |
|---|--:|---|
| HCC1937 | +0.52 | 最強 |
| HCC1395 | +0.46 | Inner TP 0.959 (n=219) |
| HCC1395_DORADO | +0.39 | 平台 replicate |
| HCC1954 | +0.34 | ⚠ 見樹：Inner TP 僅 0.43（全 TO 最低）+ polyploid 可靠性存疑 |
| H1437 | +0.23 | |
| COLO829 | +0.10 | 第 7 樣本，gap 第 2 小但同向（Inner 0.724 n=1196 / Outer 0.620 n=8380）|
| H2009 | +0.05 | ⚠ 見樹：baseline 飽和（Inner 0.93 / Outer 0.88），7/7 中最弱一票 |

[src 全表: obs18b_wilcoxon_ng2_gap_n7.json]

**Provenance fix**：原 B1 Wilcoxon 腳本從未 commit → 本 session 重寫 `obs18b_wilcoxon_ng2_gap.py`，重現 B1（W=21 p=0.0156, reproduces=True）。[src: evidence_ledger.jsonl]

---

### §2.2 支柱二（支撐）— ZAR1L/BRCA2 ASM｜SUPPORT-NULL ⭐3（單樣本 HCC1395）

**One-line Verdict**：**BRCA2/ZAR1L 位點 ASM 真實、可重現、抗 Bonferroni 且 cis-test 雙軸通過；但全基因組層級是稀有事件、無方向偏好、不能判別 TP/FP、不可一般化 → reviewer-proof supplement-null。**

**(a) 全鏈調查定案（5/29–5/31）[F]/[O]**

| 假設 | 方法 | 結果 | verdict |
|---|---|---|---|
| BRCA2 位點 ASM 真實穩健（**兩者並陳**）| 4 切法 + Bonferroni + neg control + copy_partition | **HP-axis Δβ=−0.122**（n=197, 含 allele/copy）／**純 within-somatic d_within=−0.023**（perm p=0.022，小但真）；ALLELE −0.099；chr13:32,315,128 G>A, hypo 抗 Bonferroni | ✅ 真實但**純 somatic 效應小** |
| 全基因組普遍 clear ASM？ | 嚴格多重檢定漏斗 | 39,447 TP → 30,511 unique（51,171 records）→ p<0.05 **11,830** → Bonferroni **313** → strong-ASM **172 (0.34%)** → 去重 **~26 獨立位點** | 稀有事件 |
| strong-ASM 判別 TP/FP？ | Fisher + AUC | strong-ASM 在 FP 富集 **5×**（OR=0.194, p=1.8e-28）；連續 |Δβ| AUC=**0.5049**（perm-p=0.496 中性）| ❌ NEGATIVE 且 anti-discriminative |
| somatic 驅動 read-level 聚類？（B2 subclone）| ARI vs imprinting 驗尺 | TP median ARI **0.135** vs imprinted ruler **0.758**（GNAS/RB1=1.0）；TP vs het MW p=0.922 | ❌ NEGATIVE（非特異）|
| 方向偏好？ | binomial | hypo 76 (44.2%) / hyper 96 (55.8%)，p=0.147 | ❌ 無偏好 |

[src: 20260531_ISM_aux_tag_observation_funnel_01.md / 20260529_ISM_ZAR1L_BRCA2_ASM_verification_01.md]
> **數字誠信註**：本漏斗的 `FDR 3,834 (7.5%)` 與衍生比例 `(0.61%)` 在源頭 funnel 文件 grep 不到（workflow verifier 攔截）→ **已剔除不引用**。strong-ASM 計數以 Bonferroni 172 為準。

<!-- provenance-verified: d_within −0.023 來自 research/tsg_promoter_asm_reviewer/genome_survey_v2/copy_partition_confirm.json (2026-06-03 14:06 快照; brca2_perm obs=−0.0233 p=0.022)；≈ 既有 d_drift −0.022（line 119 cis-test）。 -->
> **【2026-06-03 amendment — BRCA2 改「兩者並陳」口徑】**：tsg `copy_partition_confirm.json` 四維分解顯示 BRCA2 **HP-axis −0.122** 大半為 allele/copy（d_allele=−0.099 / d_copy=−0.11），**純 within-somatic d_within=−0.023**（perm p=0.022，小但真）。故對外引用 BRCA2 強度敘述須**溫和**：「**真實、抗 confound，但純 somatic 效應小**」，不可講「強 somatic ASM」。[src: tsg copy_partition_confirm.json @06-03 快照，project 仍 active，投稿前以定稿為準]
>
<!-- provenance-verified: d_HP −0.122/d_copy −0.11/d_within −0.023 ← copy_partition_confirm.json; chr17 within 0.142 / perm p 0.001 / BRCA2 perm p 0.024 ← survivor_permutation.json; MW p 0.6183 / signed ρ −0.083 ← fast_cnv_validation.json; HaplotagStrategy.cpp:505-516 ← longphase-s source. 全 2026-06-09 本 session 重 grep 驗證 (ledger 95)。 -->
> 🔴 **【2026-06-08 amendment + 2026-06-09 R1/R2/R3 統一口徑 — BRCA2 = subclone/copy-confounded；本 §2.2「cis-test 雙軸通過 / 真 cis-driven」口徑已過時】**（收斂稽核 wf_9e169112-573 + capstone ledger 94）：06-07 重分析把 BRCA2 (chr13:32,315,128) 列為 **predominantly subclone/copy-confounded**——**HP1-1 是 longphase-S 的 somatic SUBCLONE tag〔germline-H1 + 帶 somatic ALT, HaplotagStrategy.cpp:505-516〕，非 copy；非 CN-dosage**。HP-axis Δβ=−0.122 ≈ d_copy −0.11〔subclone vs germline 同 REF〕+ d_within −0.023〔focal cis 殘餘, perm p=0.024, **邊際·未確立**〕。816 HP-axis loci 經 Bonferroni+copy-test 後**唯一乾淨 cis = chr17:79,991,120 (TBC1D16)**（within 0.142 > HP 0.122, perm p=0.001, 單樣本）。**copy-DOSAGE 決定性非 driver**（MW p=0.6183, signed ρ=−0.083 反向 REFUTED）→ 故 **Martin-Trujillo (CN-dosage) 不 corroborate 此 reclassification**（category error；其正確角色是我們控掉的 confound 警告）。⇒ 對外**勿**再寫「BRCA2 真 cis-driven / 雙軸通過」（下文 (c) line 122 已被本 amendment supersede）；改「BRCA2 = subclone/copy 主導例（焦點 cis 單樣本不可分離），chr17/TBC1D16 = 唯一乾淨 cis exemplar」——**勿寫精確 %（split 不 robust）、勿寫「copy number」、勿寫「純 artifact」**。完整 paper-readiness + GO/NO-GO → `InterSubMod/knowledge/11_external_literature/10_paper_readiness_convergence.md`（HD-3）。⚠ tsg project 仍 active，投稿前以定稿為準。

**(b) 5/31 六-agent 紅隊 critique [O]**：核心 NEGATIVE 更穩（correction-robust：FDR 較寬鬆但 null≥TP 在兩種校正都成立）；修正 5 個正面子 claim 口徑 — **5× FP 富集中 61% 來自 chr8 hotspot**（drop-chr8 後 OR 5.16→2.84）；"FP |Δβ| larger" claim **已 RETRACTED**（p=0.137 NS）。

**(c) 6/01 BRCA2 cis-test（方法學亮點）[F]** ⚠ **【2026-06-08 SUPERSEDED — 見上 amendment：06-07 重分析後 BRCA2 = subclone/copy-confounded（HP1-1=somatic subclone tag 非 copy、非 CN-dosage；% split 不 robust），此「真 cis-driven」口徑已退役，乾淨 cis 改 chr17/TBC1D16】**：normal-anchored cis-test（normal-HP1 / tumor-HP1 / tumor-HP1-1）證 **BRCA2 真 cis-driven**（d_cis=−0.142/−0.152, d_drift≈−0.022）；高 cohesion / 高 ARI 的 germline-het 只 drift → **cohesion/ARI ≠ cis**。HP1-1 somatic silhouette 0.267 > HP1 0.119 > HP2 0.089。釐清 **MSA = 外部抽取工具 ≠ ISM binary**。[src: 20260601_ISM_BRCA2_cis_deepdive (ledger 20260601_brca2_cis_test)]

> **magnitude 演變史（口徑）**：buggy Δβ=−0.054 = MSA Level1 把 5mC+5hmC 雙列各計一次把 beta 砍半的 artifact（非 CpG-set 差異）；正確口徑 max-collapse any-mod = **−0.122**。舊 Wilcoxon p=6.09e-11 因 n 灌水 ~2.4×，**絕對值不可引用**（符號方向不受影響）。

---

### §2.3 支柱二補強 — 6/02 ASM × CN confound pilot + regime-A + LOH 診斷｜⭐3 Tier-A [F]/[O]

**One-line Verdict**：**HCC1395 單樣本下，甲基分群差異主要 mutation-linked；CN/ploidy/LOH/error/alignment confound 在所有可測通道大致排除。**

| 假說 | verdict | 關鍵數字 |
|---|---|---|
| H-A：HP-axis |Δβ| 獨立於 median CN | ✅ CONFIRM | partial ρ=**−0.055**（p=7.7e-5, n=5142）；dedup ρ=−0.039 |
| H-B：甲基無法預測 CN class | ◐ INCONCLUSIVE | dAUC=**0.046**（灰帶 0.02–0.05）；absolute AUC<0.58 實務無用 |
| H-C：HP/ALLELE 一致、excess 非 dosage | ✅ CONFIRM | MW p=0.967；HP-vs-ALLELE ρ 0.59 overall / 0.78 cnLOH |
| H-D：非 error/alignment 解釋 | ◐ CONFIRM-partial | NM/supp ρ≤0.05；⚠ **MAPQ 恆=60 是 blind spot** |
| H-E：聚類非 CN/對齊 artifact | ✅ CONFIRM | blind-ARI ≫ placebo（0.267 vs 0.015）, paired Wilcoxon **p=1.8e-43** (n=338) |
| Q1 前置 gate | ✅ PASS | β 重算 **20/20**（max abs diff=0）；PC1 ARI 0.557 / NC1 −0.005 |

[src: 20260602_ASM_CN_confound_disentanglement_pilot_01.md]

- **regime-A 乾淨子集（G1）**：75 → 62 → 23 過 length-placebo → **15 撐過完整 collider battery**（M4c 22/23=96%, M8 16/23=70%）；regime-A blind-ARI **0.229 vs het-null 0.038**（MW p=2.3e-4, Cliff δ=0.368）。⚠ 強化 collider **刷掉 SOX2/HOTTIP/SDHAF1**；**BRCA2 在 regime-A 內但未過 placebo cut（borderline, placebo 0.132>0.10）**。
- **LOH 表觀雙 haplotype 成因診斷（G2，驗證盲點假設）**：LOH-regime sig HP-axis 2170 位點（cnLOH 816 / gainLOH 1343 / lossLOH 11）；分層抽樣 79 = **57 (72%) self-phasing artifact / 14 (18%) candidate subclone / 7 (9%) CN-regression / 1 (1%)**。診斷非強制歸類（強制歸類=V5 Layer 1.5 已失敗）。
- **口徑校正**：分析 TP set = 30,490 ClairS-TO calls（recall 77.3%），**非** SEQC2 truth 39,447。master 56,320 records / unique 34,155（TP 30,511 + FP 3,644）。

[src: 20260602_ASM_cluster_characterization_regimeA_LOH_diagnostic_01.md]

**天花板**：全部單樣本 ⭐3 Tier-A，**M9 cross-sample 結構性不可達**（僅 HCC1395 有 ground-truth）→ 不可宣稱跨樣本一般性。

---

### §2.4 支柱三（methods）— LOSO sample-level NEGATIVE｜⭐⭐⭐⭐ L2 [F]

本週無新分析，但 T4 handoff 紀律強化其角色：HCC1395 in-dist +0.02236 vs **LOSO held-out −0.00012**（drop +0.02248 = **100% sample-level circularity**）；5-sample mean −0.00004；Wilcoxon p=0.125, 0/5 pos 4/5 neg → **DIRECTION_NEGATIVE**；HCC1954 transfer −0.377。論文 methods 段的誠實方法學貢獻。[src: CURRENT_FOCUS.md]

---

## §3 本週其他活動（非論文三軸）

### §3.1 外部 handoff / methods 設計（T4，全 in_progress）[O]/[U]
- **HP3/unphase 甲基 phasing-rescue** = **DESIGN-ONLY**，Phase A GO/NO-GO gate 未跑（MD1-MD10 實測值欄全 `—`）→ 「甲基能救 phasing」**目前無實測支持**。基線（既有行為）：120,834 somatic 位點 read-instance；H3 救 42%；untag 1,777,511 不變；H0 65,033（read-instance, Tmode build）。
- **methyl-assisted untag 二次分類** = feasibility-landscape ⭐2 L4（純文獻 + 既有實作盤點，無實測）。
- **BRCA2/ZAR1L TSG promoter ASM reviewer-response** = POSITIVE 但**單樣本 n=1**（germline-HP1 β=0.215 n=63 vs somatic-HP1-1 β=0.093 n=47, Δβ=−0.122, Wilcoxon p=6.1e-11, Cohen d=−0.394）。[src: 20260527_HCC1395_TSG_promoter_ASM_reviewer_response_01.md]
- **LongPhase-S（HKU 黃耀廷 lab）第三方 audit** = ⭐⭐ **2/5 Major revision**（3 fatal gap；R1/R2/R3 = ⭐2/2/2.5；completion 35-40%）。[src: 20260527_LongPhase-S_revision_response_3rd_party_audit_01.md]

### §3.2 治理 / infra（T5，間接護航論文）[F]
- **disk cleanup 回收 big7 ≈623 GiB**（278.42+278.42+112.23 GB = 669.07 GB；df 97%→96%，1.4T→2.0T free）；只刪 3 個可重生 TO BAM，**全保留 LOH.bed/vcf provenance**；canonical 14 tagged BAM **3.37 TiB KEEP**（ASM live working set）；big8 唯讀不可動。[src: FINAL_CLEANUP_LEDGER.md]
- **HTML 捏造事件 postmortem**：同 session 30 分鐘內**純文字規則失效 3 次、方向皆與真值相反**（捏造 BRCA2=0.572 真值 0.866；null median 0.578 真值 0.974；V4 75% 真值 95%）→ 落地 **§13 三層機械防線**（物理隔離分析與寫報告 batch）。[src: 20260601_fabricated_metric postmortem]
- **11-agent harness audit** 裁決 **0 個新外部工具**（19 條去重 ADD6/MODIFY5/ALREADY3/SKIP_HYPE2/SKIP_RISK3）；最高破口 = `pre_tier_upgrade_check` 未 wired。[src: harness_audit_dashboard_design doc]
- git log：今日 commit `7a4cc6f`/`37e1511` 採用 obra/superpowers 流程紀律借鑑 B1-B8（源頭在 CLAUDE.md，非 T5 三檔）。

---

## §16 下一步行動（Layer 4）

- ⭐⭐⭐ **Decisive**：啟動 **G6 framework draft** 寫作（phasing 主軸 B+ 證據定稿）；full-Grade-A 負控（~25-50hr C++ flag-on）**延後至 reviewer 要求**（本週決議）。
- ⭐⭐ Operational：ASM 支撐軸以收斂口徑（real / non-directional / non-discriminative）寫進 draft 的 characterization + supplement-null 段。
- ⭐⭐ Operational：CURRENT_FOCUS line 48 仍存舊「方向 POSITIVE」表述（line 25 已收窄）→ 寫作前統一為收斂口徑。
- ⭐ Maintenance：stale queue H015-018（filter 已 DEAD）待 /pivot-direction 降權；`pre_commit_compile_check`/`kb_schema_check` neutering + stale marker 待 Hard Gate 決策。

---

## §17 教授可能提問 + 預備回答（Layer 4）

**⭐⭐⭐ 必問（Must-Answer）**
1. *「n=7 是 7 條獨立樣本嗎？」* → 6 條 cell line + 1 平台 replicate（HCC1395 DORADO）；方向全一致，最弱兩票 H2009（baseline 飽和）+ COLO829（跨平台）仍同向。
2. *「為什麼還不是 full Grade A？」* → 自參考軸 circularity（Inner 用 HPFineN_HP1S=HP1-1 count）完整負控需 ~25-50hr C++ flag-on 全 7-sample 重跑；機制已有 X3+V6C 4-sample chr19 支持，可先寫作、reviewer 要求再補。
3. *「ASM 既然不能 filter，留它幹嘛？」* → 當 characterization + supplement-null；BRCA2 是真 cis 案例（normal-anchored cis-test 證實），與全域稀有並存——正是誠實的 read-level epigenetic context 故事。

**⭐⭐ 可能問（May-Ask）**
- *「ASM 5× FP 富集不就是有判別力？」* → 那是 LOH→single-hap→low-cov→極端 baseline 的回歸極值 artifact（ONE regime），且 61% 來自 chr8 hotspot；連續 |Δβ| AUC=0.505 中性才是判別力的正確讀數。
- *「6/02 confound pilot 排除了所有 confound 嗎？」* → 可測通道 H-A/C/D/E 排除、H-B 灰帶；唯一 blind spot 是 MAPQ 恆=60（無變異無法測 low-MAPQ confound）。單樣本天花板，跨樣本結構性不可達。

---

## §18 敘述腳本（如何對 PI 解釋）

**60 秒版**：
> 「這篇論文有三根支柱。**主軸 phasing 本週往前一大步**——加進 COLO829 後跨 7 個樣本方向全一致，統計強度到了 Grade A 門檻。**支撐軸 ASM 我這週做的不是找更多正向，而是把它的弱點全部自己先證偽**：證明它不能判別真假突變、排除了它的 magnitude 是 copy-number dosage 驅動（個別最強位點如 BRCA2 仍被 copy-partition 混淆，焦點 cis 單樣本不可分離）、也排除比對假象——所以它在論文裡是乾淨的 characterization 加一個 reviewer 打不動的 null。**第三根（LOSO 負面結果）**本來就是方法學賣點。整體上，這是把論文從 B+ 推向 A、同時把所有會被 reviewer 質疑的點提前關掉的一週。」

**核心 framing**：本週價值 = **收斂 + 釘死邊界**，不是擴張。把「ASM 很誘人但會被打」的 over-claim 自己先拆掉，反而讓論文更硬。

---

## 數據可信度與來源（provenance footer）

- 全 metric 數字經 2026-06-02 gather→verify workflow（`wf_02d8afb0-e13`, 12 agents, 596s）雙階段對抗 grep 驗證：T1 19/19✓ / T2 27/28（剔除 1 grep-fail）/ T3 32/32✓ / T4 21/21✓ / T5 16/16✓ / T6 11/11✓。
- 口徑提醒：read-instance vs unique、Tmode vs paired、FDR vs Bonferroni 不可並列；ASM 全為 HCC1395 單樣本 ⭐3 Tier-A（⭐4 需 COLO829，標 off-goal）。
- git HEAD（撰寫時）：37e1511。報告 provenance 鏈見各 data_sources 檔。
