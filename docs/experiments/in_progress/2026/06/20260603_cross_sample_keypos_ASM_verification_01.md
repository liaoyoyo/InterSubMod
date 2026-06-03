---
title: 跨樣本關鍵位點 ASM 驗證 — HCC1395 38 位點 × 6 癌症樣本
date: 2026-06-03
anchor_commit: 12b654a
branch: feat/cis-asm-pipeline
task_type: A pilot (targeted cross-sample + genome-wide expansion)
tier: ⭐3 (bordering ⭐4)
partial_flag: true
partial_scope: "6 樣本 × 38 targeted key positions + 6 樣本 × 2000 random TP somatic SNVs (genome-wide subset, 非全 somatic)"
data_sources: research/tsg_promoter_asm_reviewer/genome_survey_v2/cn_confound/cross_sample/cross_sample_synthesis.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/cn_confound/cross_sample/HCC1395_gwasm.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/cn_confound/cross_sample/COLO829_gwasm.json
deliverable: docs/experiments/in_progress/2026/06/20260603_cross_sample_keypos_ASM_verification_01.standalone.html
evaluator_verdicts: "targeted separation-correctness:PASS; targeted over-claim:NEEDS_WORK→fixed; genome-wide reproducibility:PASS+1fix"
ledger: 20260603_cross_sample_keypos_ASM_verification, 20260603_genomewide_asm_rate_crosssample
---

# 跨樣本關鍵位點 ASM 驗證

> **互動可視化交付物**：`InterSubMod/docs/experiments/in_progress/2026/06/20260603_cross_sample_keypos_ASM_verification_01.standalone.html`（212 KB，3 圖 7 表 + genome-wide 章節，數字 §13 layer-A 由 JSON 注入）。本 .md 為 provenance 錨定記錄；所有 metric 可在 `data_sources` 的 JSON grep 到。

## TL;DR

1. **同位點 somatic 復發 = 0/38** —— 其他 5 個癌症樣本在 HCC1395 的 38 個關鍵 ASM 位點全部 **0 個** somatic call（各癌 private mutation）→ 證實 ASM 的 **somatic 因果性**（若是 CN/技術假象應跨樣本復發）。
2. **somatic ASM 現象跨癌種復現** —— 各樣本用自己的 private somatic 突變，全 **6/6 樣本 excess-over-null > 0**（乳腺/肺/黑色素瘤皆正），但這是**單一管線的「現象複製」非獨立管線驗證** → tier 封頂 ⭐3 bordering ⭐4。

## 1. 場景與資料

| 樣本 | 癌種 | tagged BAM | methylation | somatic VCF |
|------|------|:---:|:---:|:---:|
| HCC1395 | 乳腺（SEQC2 ref） | 260G | MM/ML | ✓ |
| HCC1937 | 乳腺（BRCA1-mut） | 417G | ✓ | ✓ |
| HCC1954 | 乳腺 | 225G | ✓ | ✓ |
| H1437 | 肺腺癌 | 215G | ✓ | ✓ |
| H2009 | 肺腺癌 | 290G | ✓ | ✓ |
| COLO829 | 黑色素瘤 | 89G | ✓ | ✓ |

- **HCC 系列同 Hamon Cancer Center 但不同病人**（非同細胞株）。
- BAM 用 **`HP:Z:` 字串 tag**（值 `1`/`2`/`1-1` somatic 子單倍型），非標準 `HP:i:`。
- 38 關鍵位點 = 37 credible_loci（HP-axis somatic-controlled）+ flagship **BRCA2/ZAR1L chr13:32315128**（ZAR1L = BRCA2 同位點，ZAR1L 為 nearest gene）。

## 2. 方法

- 複用已驗證 pysam 抽取（54 modkit crossval Pearson=1.0）：per-read 5mC β（mean over CpGs in ±600bp），依 `HP:Z` tag 分群；somatic-controlled HP-axis Δβ = β(subhap) − β(main)。
- **Targeted**（script 57 + 確定性合成 58 v2）：38 位點 × 6 樣本；雙軸（somatic / germline）+ sign-concordance。
- **Genome-wide**（script 59）：各樣本自己 `filtered_snv_tp.vcf.gz` 隨機 N=2000；HP-label-shuffle null（20 perm）估噪音底；**excess-over-null = 真訊號**。
- **HTML**（script 60）：§13 layer-A，數字從 JSON 注入。**Ledger** append（script 61）。

## 3. Targeted 結果（38 位點 × 6 樣本）

- 同位點 somatic 復發 **0/38**；HCC1395 somatic 重疊 33/38（28 somatic-controlled ASM），其他 5 樣本全 0。
- **28 HCC1395-private somatic ASM**（含 BRCA2/ZAR1L HP-axis Δβ≈−0.19）。
- BRCA2 逐樣本：只有 HCC1395 有真 somatic 子單倍型 + ASM；COLO829/HCC1954 = germline het 雙倍型、H1437/HCC1937 = LOH 單倍型、H2009 = 4-read 假性 subhap（Δ≈0）。
- germline 軸復發 8 個 → 嚴格方向檢查拆：**4 sign-concordant 候選 imprinting**（HERC6/HOTTIP/LOC101927914/LOC124903622）+ **4 方向相反 NON-imprinting**（SOX2/LINC00689/LOC283028/DUX4L33）。

## 4. Genome-wide 復現結果（各樣本自己 private somatic，N=2000）

| 樣本 | 癌種 | excess-over-null |
|------|------|:---:|
| HCC1395 | 乳腺 | +0.241 |
| HCC1954 | 乳腺 | +0.196 |
| HCC1937 | 乳腺 | +0.171 |
| H2009 | 肺 | +0.151 |
| H1437 | 肺 | +0.150 |
| COLO829 | 黑色素瘤 | +0.101 |

- **全 6/6 excess-over-null > 0**（mean 0.168，CV 0.26，3 癌種皆正）→ somatic ASM 現象跨癌種復現。
- **必看 excess 不可用 raw rate**：median|Δβ| 全 <0.10（固定門檻 raw rate 受小 N 噪音灌水）；null 隨覆蓋變動（COLO829 null 最高因 n_evaluable=802 最低）。
- 方向多輕微 hyper（52–58%），HCC1954 唯一輕微 hypo。imprinted DMR 正控（EXPLORATORY）PEG3 強（HCC1395 +0.72 / HCC1937 −0.93，反映各自 LOH 保留親本 allele 不同）。

## 5. 驗證（generator–evaluator 分離）

- **Targeted evaluator 1（分離正確性）：PASS** — 標籤全 data-grounded，交叉驗證 6 cells 0 transcription error。
- **Targeted evaluator 2（over-claim）：NEEDS_WORK → 已修** — 抓到「imprinting」誤用（4 hit 中 3 個跨樣本方向相反）→ 改方向不可知命名 + sign-concordance。
- **Genome-wide evaluator（復現 + null honesty）：PASS + 1 必補** — 已補「單一管線 caveat」入 HTML。

## 6. 限制（誠實標註）

1. SEQC2 truth 只有 HCC1395（其他 5 樣本 TP VCF 是 caller 輸出，無正交 truth）。
2. targeted「不復發」部分受 LOH/低 N 使 germline 軸不可測所限（非純生物缺席）；IGF2 cnLOH 使已知 imprinting 被遮蔽。
3. **6 樣本共用同一 ClairS-TO/longphase caller + 同一 HP-axis 估計法 → 現象複製非獨立管線驗證，共用系統性偏差未排除 → tier 封頂 ⭐3**。
4. genome-wide rate 用 window-mean 5mC，與 targeted credible-loci 的 paired-CpG MAX-collapse Wilcoxon 口徑不同，方向可比、magnitude 不可並列。
5. genome-wide 用各樣本 2000 random TP somatic（非全 somatic）；未出 per-sample permutation p（僅 point-estimate excess）。

## 7. 後續（⭐4 升級路徑）

- 獨立管線（非 ClairS-TO）或正交 truth 複製 + per-sample permutation p on excess。
- COLO829 melanoma 已納入 genome-wide ASM（原 ⭐4 blocker 之一首次有數據）。
- 作為 phasing 論文中 ASM 真實性的 methods 支撐（somatic-driven + 位點 private + 現象跨癌種復現）。

## See Also

- Memory: `project_cross_sample_asm_reproducibility`（+ [[project_asm_cn_confound_pilot]] / [[project_zar1l_brca2_asm_verification]] / [[feedback_asm_allele_axis_baseline_confound]]）
- Scripts: `research/tsg_promoter_asm_reviewer/scripts/57-61`
- Data: `research/tsg_promoter_asm_reviewer/genome_survey_v2/cn_confound/cross_sample/`
