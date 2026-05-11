---
id: ism-kb-09-conclusions-concluded-negative
name: "Concluded NEGATIVE / NO-GO 目錄"
description: "已證偽的 20+ 方向清單；別重新調查。含 O9-O13、Option C、Fine-Pairwise、TO Germline FP、Beyond-AUC、Wave 3 LOH 等。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: reference
verified_scope: "NEGATIVE conclusions against MEMORY Concluded section"
related_ids:
  - ism-kb-09-conclusions-index
  - ism-kb-09-conclusions-characterization-only
  - ism-kb-02-samples-hcc1954
  - ism-kb-03-pipelines-tumor-only
  - ism-kb-07-derived-features-fisher-frac-sig
  - ism-kb-09-conclusions-negative-01-readparser-germline-hp-only
  - ism-kb-09-conclusions-negative-02-r1-global-sample-asm
  - ism-kb-09-conclusions-negative-03-clairs-to-verdict
  - ism-kb-09-conclusions-negative-04-z3-internal-exploration
  - ism-kb-09-conclusions-negative-05-o11-heterogeneity
  - ism-kb-09-conclusions-negative-06-o12-loh-methylation-scenarios
  - ism-kb-09-conclusions-negative-07-o13-cross-region
  - ism-kb-09-conclusions-negative-08-germline-fp-nogo
  - ism-kb-09-conclusions-negative-09-o9-fn-characterization
  - ism-kb-09-conclusions-negative-10-to-pure-independent
  - ism-kb-09-conclusions-negative-11-fine-pairwise
  - ism-kb-09-conclusions-negative-12-beyond-auc-exhaustion
  - ism-kb-09-conclusions-negative-13-option-c-dual-path
  - ism-kb-09-conclusions-negative-14-wave3-loh-stratification
  - ism-kb-09-conclusions-negative-15-feature-design-r1r5
  - ism-kb-09-conclusions-negative-16-qs-redesign-evidence
  - ism-kb-09-conclusions-negative-17-to-feature-deep-study
  - ism-kb-09-conclusions-negative-18-ism-af-discriminability
  - ism-kb-09-conclusions-negative-19-to-fp-provenance
  - ism-kb-09-conclusions-negative-20-paired-f1-filter
tags: [conclusions, negative, no-go, concluded, index]
canonical_paths: [09_conclusions/03_concluded_negative.md]
alias_paths: []
---

# Concluded NEGATIVE / NO-GO 目錄

- 一句結論：20+ 方向已證偽，**別重新調查**；含失敗原因與替代方向
- 適用對象：AI agent（避免重複）、研究決策前檢核
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  grep "^- \[" /bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory/MEMORY.md | \
    grep -i "negative\|no-go\|concluded"
  ```

---

## 🔍 關鍵字快速索引（2026-04-23 新增）

| 關鍵字 | 對應 N 編號 + MEMORY ID |
|--------|-------------------------|
| HPFineNGroups flag=on | N1（`readparser_germline_hp_only_phase1_negative`）|
| Sample ASM | N2（`r1_global_sample_asm_negative`）|
| ClairS-TO Verdict | N3（`clairs_to_verdict_pilot`）|
| Z3 / HCC1954 amplicon | N4（`z3_internal_exploration_negative`）|
| Epipolymorphism / heterogeneity | N5（`O11_heterogeneity_negative`）|
| L2 collider bias / LOH methylation | N6（`O12_loh_methylation_scenarios`）|
| Cross-region methylation | N7（`O13_cross_region_negative`）|
| TO Germline FP | N8（`germline_fp_identification_nogo`）|
| FN rescue HP-free | N9（`O9_fn_characterization_nogo`）|
| TO-pure | N10（`to_pure_independent_modeling_negative`）|
| Fine-Pairwise 距離 | N11（`fine_pairwise_negative`）|
| Beyond-AUC / 特徵空間耗盡 | N12（`beyond_auc_exhaustion_confirmed`）|
| Option C / HP-free dual path | N13（`option_c_dual_path_negative`）|
| Wave 3 LOH / cnLOH | N14（`wave3_loh_stratification_closure`）|
| CramersV 2×2 缺陷 | N15（`feature_design_limitations_r1r5`）|
| QS TO 隨機 | N16（`qs_redesign_evidence`）|
| TO RF ceiling | N17（`to_feature_deep_study`）|
| ISM vs AF 全量平坦 | N18（`ism_af_discriminability`）|
| TO FP provenance | N19（`to_fp_provenance`）|
| Fisher_Frac_Sig / F1-filter | N20（`paired_f1_filter_abandoned`）|

**用法**：查某主題是否已證偽 → 在此表 grep 關鍵字 → 跳對應 N 編號 → 讀詳情。

---

## NEGATIVE 清單（v0.5 ADR 化；點 N 編號讀完整 Context/Decision/Method/Result/Consequences）

每項格式：**[N 編號] 標題 — 一句結論 — [完整 ADR 連結]**

### N1 — [ReadParser `--germline-hp-only` Phase 1 CONDITIONAL NEGATIVE](negative/N01_readparser_germline_hp_only.md)
- 2026-04-21｜機制正確但 HCC1395 TO 全量無 TSV 特徵 AUC ≥0.02 改善；**副作用**：HPFineNGroups N≥3 在 flag=on 下消失

### N2 — [R1-Global Sample ASM NEGATIVE](negative/N02_r1_global_sample_asm.md)
- 2026-04-21｜SampleASM_Delta residualized AUC 0.527 [CI 0.520, 0.533]；CL-025a ⭐3→⭐2

### N3 — [ClairS-TO Verdict Pilot NEGATIVE on F1](negative/N03_clairs_to_verdict.md)
- 2026-04-20｜內部校準 POSITIVE 但 Verdict_Germline 100% 落 LowQual；主升級路徑改 Wakhan/SAVANA

### N4 — [Z3 Internal Exploration NEGATIVE](negative/N04_z3_internal_exploration.md)
- Z3 內 12 特徵 AUC<0.61；HCC1954 amplicon 已知

### N5 — [O11 Heterogeneity NEGATIVE](negative/N05_O11_heterogeneity.md)
- epipolymorphism AUC 全因 n_reads confound

### N6 — [O12 LOH 甲基化場景 NEGATIVE](negative/N06_O12_loh_methylation_scenarios.md)
- AlleleDelta=AF confound；**L2 collider bias**（方法論教訓：殘差化必走 within-group OLS + AF-bin 交叉）

### N7 — [O13 跨區域甲基化 NEGATIVE](negative/N07_O13_cross_region.md)
- shared read count confound；分層後消失

### N8 — [TO Germline FP NO-GO](negative/N08_germline_fp_nogo.md)
- G1-G7 全 AUC<0.64；FP removal=0%

### N9 — [O9 FN NO-GO](negative/N09_O9_fn_characterization.md)
- HP-free AUC<0.53；FN≡TP in methylation space

### N10 — [TO-pure NEGATIVE](negative/N10_to_pure_independent.md)
- caller_af=0.654 超越全 ISM

### N11 — [Fine-Pairwise NEGATIVE](negative/N11_fine_pairwise.md)
- 6 pairwise 距離全無效；特徵空間耗盡

### N12 — [Beyond-AUC 耗盡確認](negative/N12_beyond_auc_exhaustion.md)
- pure methylation 全≤0.58；HPFineNGroups somatic marker 確認

### N13 — [Option C 雙路 NEGATIVE](negative/N13_option_c_dual_path.md)
- HP-free AUC=0.564；所有信號來自 HP tags

### N14 — [Wave 3 LOH 分層 CLOSED](negative/N14_wave3_loh_stratification.md)
- Non-LOH 無突破；cnLOH Simpson's Paradox

### N15 — [Feature Design R1-R5](negative/N15_feature_design_r1r5.md)
- CramersV 93%零=2×2 缺陷；**identifiability 根因**

### N16 — [QS TO 失效](negative/N16_qs_redesign_evidence.md)
- QS AUC=0.497 隨機；需根本重設計

### N17 — [TO Feature Deep Study](negative/N17_to_feature_deep_study.md)
- RF ceiling 0.69-0.77；甲基化全無效

### N18 — [ISM Methylation vs AF](negative/N18_ism_af_discriminability.md)
- 全量效果平坦

### N19 — [TO FP Provenance](negative/N19_to_fp_provenance.md)
- TO FP 無法 ISM 過濾；FN rescue 方向後被 N9 否決（**完整封閉圈**）

### N20 — [Paired F1-filter 方向放棄](negative/N20_paired_f1_filter.md)
- 2026-04-21｜Fisher_Frac_Sig CI 跨隨機；characterization-only

---

## 核心 Pattern：為何 NEGATIVE？

大多 NEGATIVE 源自以下 pattern：

### Pattern 1：Confound with AF or coverage
- O11, O12, O13, Heterogeneity → 信號被 nuisance 解釋

### Pattern 2：Collider bias（residualization on AF）
- L2 residualize, O12 → 產生假信號

### Pattern 3：FP removal=0%
- TO Germline FP, Verdict → 有判別力但無實務過濾力

### Pattern 4：HP tag artifact
- HPFineNGroups flag=on 下消失、Option C HP-free AUC 降
- → 許多「信號」來自 LongPhase HP tag 人工分組

### Pattern 5：特徵空間耗盡
- Fine-Pairwise, Beyond-AUC → pure methylation 全≤0.58

---

## 對 AI agent 的建議

看到以下關鍵字時，**先查本目錄**：
- "能不能做 XX filter"
- "有 pattern 似乎能當 signal"
- "AUC > 0.55 的特徵"
- "新發現的 characterization"

**大機率已經試過**；除非使用者明確說要挑戰既有結論。

---

## 相關

- Characterization（類似但 positive）：[02_characterization_only.md](02_characterization_only.md)
- MEMORY 索引：`/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory/MEMORY.md`
- Landscape 總索引：[04_research_landscape_index.md](04_research_landscape_index.md)
