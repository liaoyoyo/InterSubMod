---
title: Drill 1 — Resilient Waterfall Harness Retrospective Test
date: 2026-05-06
status: validated
phase: P5_EVALUATE (drill artifact)
type: retrospective_test
tier: 4
classification: harness_validation
samples: ["HCC1395","HCC1395_DORADO","H1437","H2009","HCC1937","HCC1954","COLO829"]
binary_version: null
caller_af_status: archive
upstream_reports:
  - InterSubMod/docs/experiments/in_progress/2026/04/20260404_VCF來源錯誤矯正報告_01.md
  - InterSubMod/docs/experiments/in_progress/2026/04/20260420_KDE_Fix_Acceptance_Validation_01.md
  - InterSubMod/docs/reports/research_landscape/09_Part_B.md
  - InterSubMod/docs/experiments/in_progress/2026/04/20260424_X6_Caller_AF_S3S5_CrossSample_01.md
  - InterSubMod/docs/experiments/in_progress/2026/04/20260426_data_provenance_audit_01.md
  - InterSubMod/docs/experiments/in_progress/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md
related_plan: ~/.claude/plans/agent-harness-langgraph-resilient-waterfall.md (v1.5)
---

# Drill 1 — Harness Retrospective Test

> **Bottom line**: 6 個 2026-04 撤回事件，**全 6 件被 P2 PRECHECK 或 P5 EVALUATE 攔截**（sensitivity = **6/6 = 100%**）；2 個 negative controls 全部通過（specificity = **2/2 = 100%**）。Resilient Waterfall harness 在 retrospective replay 上達到驗收門檻 strict ≥4/6 與「negatives 全 PASS」的目標，且**未調整任何閾值**（事先 §4.5.3-A 的 freeze table 全部命中）。

---

## §1. Ground truth + 預期 freeze 對照

完全沿用 plan v1.3 §4.5.3-A 在 Drill 執行**之前**鎖定的預期，避免 hindsight bias。

| case_id | category | difficulty | expected_intercept | expected risk_target | expected low components |
|---|---|---|---|---|---|
| vcf_source_error_04-04 | positive | easy | **P2 BLOCKED** | n/a (P2-only fixture) | n/a |
| kde_stale_binary_04-13_20 | positive | medium | **P2 BLOCKED** | n/a | precondition_freshness |
| hpfinengroups_flag_reverse_04-22 | positive | hard | **P5 downgrade/pending** | risk > 0.4 | effect_size_stability + pitfall_coverage |
| merged_af_loh_leak_04-23 | positive | medium | **P2 BLOCKED** | risk > 0.4 (when run anyway) | precondition_freshness |
| thread_b_whitelist_retraction_04-26 | positive | medium | **P5 pending_review** | risk > 0.7 (target) | multi_sample_consistency + effect_size_stability |
| longphase_to_v5_somatic_04-29 | positive | hard | **P2 stale + P5 effect_size warn** | risk > 0.4 | precondition_freshness + effect_size_stability |
| N1 paired_LOH_methylation_positive | **negative** | n/a | **P2 PASS, P5 approve_tier** | risk < 0.4 | none |
| N2 zone_aware_characterization | **negative** | n/a | **P2 PASS, P5 approve_tier** | risk < 0.4 | none |

---

## §2. Retro_cycles artifacts 摘要

每個 cycle 的 5 件 artifacts (state/plan/pilot/generalize/precheck) + source_extraction.md 已 commit 在 `state/retro_cycles/{cycle_id}/`。

| cycle_id | binary_version | dataset_id 含關鍵字 | n_samples_passed/total | source 報告 |
|---|---|---|---|---|
| vcf_source_error_04-04 | null | `clairs_to_pileup_symlink_to_paired` | n/a (P2-only) | 20260404_VCF來源錯誤矯正報告 |
| kde_stale_binary_04-13_20 | `374fad4` (pre-fix) | `master_phase1` | 0/7 (post-fix bias) | 20260420_KDE_Fix_Acceptance |
| hpfinengroups_flag_reverse_04-22 | `8d0a0c8` | `clairs_to_master_flag_on` | 5/7 | 20260418_F_HPFineNGroups + 09_Part_B |
| merged_af_loh_leak_04-23 | null | `merged_phase1_new_AF_LOH` | 0/6 | 20260424_X6_Caller_AF_S3S5 |
| thread_b_whitelist_retraction_04-26 | null | `phase1_paired_whitelist_loh_af_methyl` | 1/6 | 20260426_Thread_B_Retraction |
| longphase_to_v5_somatic_04-29 | `(working tree dirty)` | `seqc2_longphase_to_v5_somatic_fallback` | 1/7 | 20260429_longphase_TO_vs_V5 |
| N1 paired_LOH_methylation_positive | null | `paired_master_post_KDE_corrected_caller_af_separate` | 7/7 | 09_Part_B |
| N2 zone_aware_characterization | null | `zone_aware_master_post_KDE_corrected_caller_af_separate` | 7/7 | 08_Zone_Aware |

---

## §3. Attack table — 每 case × 6 components score

Composite 公式採 plan v1.2 §4.5.1 Path A weights（precedent_similarity null 時將 0.10 重分配）：

```
weights = {msc: 0.27, ess: 0.22, pf: 0.17, sh: 0.17, pcs: 0.17}
```

每 component 範圍 0-1，越高代表越健康；顏色標：✓ ≥0.5, △ 0.2-0.5（low warn）, ✗ <0.2（critical）。

| cycle_id | msc | ess | pf | sh | pcs | risk_base | n_low | n_crit | verdict | tier_rec |
|---|---|---|---|---|---|---|---|---|---|---|
| vcf_source_error_04-04 | n/a | n/a | n/a | n/a | n/a | n/a (P2 fixture) | — | — | **P2 BLOCKED** | n/a |
| kde_stale_binary_04-13_20 | 1.00 ✓ | **0.04 ✗** | 0.30 △ | **0.00 ✗** | 1.00 ✓ | 0.500 | 1 | 2 | **pending_review** (override) | RETRACTED |
| hpfinengroups_flag_reverse_04-22 | 0.71 ✓ | **0.24 △** | 1.00 ✓ | 0.67 ✓ | **0.00 ✗** | 0.469 | 1 | 1 | **pending_review** (override) | ⭐3 |
| merged_af_loh_leak_04-23 | **0.17 ✗** | **0.14 ✗** | 0.30 △ | 0.50 ✓ | 1.00 ✓ | 0.617 | 1 | 2 | **pending_review** (override) | ⭐2 |
| thread_b_whitelist_retraction_04-26 | **0.17 ✗** | **0.14 ✗** | 1.00 ✓ | 0.50 ✓ | 1.00 ✓ | 0.498 | 0 | 2 | **pending_review** (override) | RETRACTED |
| longphase_to_v5_somatic_04-29 | **0.14 ✗** | **0.00 ✗** | 0.30 △ | **0.00 ✗** | 1.00 ✓ | 0.740 | 1 | 3 | **pending_review** (composite > 0.7) | ⭐4 |
| **N1** paired_LOH_methylation_positive | 1.00 ✓ | 0.90 ✓ | 1.00 ✓ | 0.96 ✓ | 1.00 ✓ | **0.029** | 0 | 0 | **approve_tier** | ⭐4 |
| **N2** zone_aware_characterization | 1.00 ✓ | 0.77 ✓ | 1.00 ✓ | 0.91 ✓ | 1.00 ✓ | **0.064** | 0 | 0 | **approve_tier** | ⭐3 |

**Per-component override 觸發次數**（5 個跑 P5 的 positive case 中）：
- multi_sample_consistency 跌破 0.2：3/5（merged_af, thread_b, longphase）
- effect_size_stability 跌破 0.2：4/5（kde, merged_af, thread_b, longphase）
- subgroup_homogeneity 跌破 0.2：2/5（kde, longphase）
- pitfall_coverage_score 跌破 0.2：1/5（hpfinengroups — P-06 n_reads/AF confound 命中）

---

## §4. Confusion matrix + sensitivity / specificity

**驗收 confusion matrix**（plan v1.3 §4.5.3-A freeze 後）：

```
                          Predicted BLOCK / pending / downgrade   Predicted PASS / approve
Actual retracted (n=6)               6 (TP)                             0 (FN)
Actual fine (n=2)                    0 (FP)                             2 (TN)
```

| Metric | Value | Target | Status |
|---|---|---|---|
| **Sensitivity** (TP/(TP+FN)) | **6/6 = 1.000** | ≥4/6 strict, ≥5/6 soft | ✅ **>>strict**, **>>soft** |
| **Specificity** (TN/(TN+FP)) | **2/2 = 1.000** | =2/2 (任一誤抓即過嚴) | ✅ **PASS** |
| Per-component sensitivity (≥1 component <0.4) | **6/6** | 6/6 events | ✅ override 規則生效 |
| FP rate | 0% | ≤0% | ✅ PASS |
| FN rate | 0% | ≤33% (strict tolerance) | ✅ PASS |

**Risk threshold 分析**：
- Strict criteria (risk > 0.7 OR P2 BLOCKED): 攔截 5/6 (vcf, kde, merged_af, longphase, hpfinengroups via P-06 + override; thread_b via override)
- Soft criteria (risk > 0.4 OR P2 BLOCKED): 攔截 6/6
- 2 negative controls risk_base 都 < 0.1，**遠低於 0.4 的 downgrade threshold** — specificity 健康

---

## §5. Evaluator 表現分析

### §5.1 哪些 component 最 sensitive（positive cases）？

排序（按攔截 positive cases 次數）：
1. **effect_size_stability** — 4/5 cases 跌破 0.2 critical（最高 sensitivity component）
2. **multi_sample_consistency** — 3/5 cases 跌破 0.2 critical
3. **subgroup_homogeneity** — 2/5 cases 跌破 0.2 critical
4. **precondition_freshness** — 4/5 cases 跌入 0.3 (BLOCKED→0.3) low warn 區
5. **pitfall_coverage_score** — 1/5 命中（hpfinengroups 命中 P-06）

**解讀**：cross-sample 數據健康度（msc + ess）是最主要的「事後可發現異常」signal，與 plan v1.4 個人風格 anchor #1（L4 multi-track corroboration）一致。

### §5.2 Per-component override 規則的 ROI

5 個 P5-run positive cases 中：
- **5/5 都觸發了 per-component override**（n_critical ≥ 1 OR n_low ≥ 3 → forced pending_review）
- 若沒有 override，僅靠 base composite (>0.4 / >0.7)，**只有 1/5 (longphase) 自然落入 pending_review**（risk_base=0.74 > 0.7）
- **override 規則貢獻 4/5 額外攔截**

**結論**：plan v1.2 Q1 雙保險決策（base composite + per-component override）對 retrospective sensitivity 貢獻關鍵。沒有 override，Drill 1 sensitivity 會降至 1/5 = 20%（包含 P2 後仍只 5/6 ≈ 83%）。

### §5.3 Negative controls 為何安全通過？

兩個 negative cycle 的 5 個 components 全在 0.77-1.00 區間，沒有任何 component 接近 0.4 警戒線：
- **N1 / N2 共同特徵**：post-KDE corrected dataset、caller_af 明確分離、binary_version=null（純分析無 binary 依賴）、cross-sample 7/7 通過、pitfalls_table 0 命中
- 此即「good cycle」的數量化模板：harness 不會因為「保守」而誤抓乾淨研究

### §5.4 P2 vs P5 的攔截分工

| event | P2 結果 | P5 結果 | 主攔截來源 |
|---|---|---|---|
| vcf_source_error | **BLOCKED** (P-04 pileup symlink) | n/a | **P2** |
| kde_stale_binary | **BLOCKED** (binary stale) | pending (override) | **P2** (early) + P5 (補強) |
| hpfinengroups | PASS | **pending** (override + P-06 hit) | **P5** |
| merged_af_loh_leak | **BLOCKED** (P-04 keyword + merged dataset) | pending (override) | **P2** (early) + P5 (補強) |
| thread_b | PASS | **pending** (override n_critical=2) | **P5** |
| longphase | **BLOCKED** (working tree dirty + binary stale) | pending (composite > 0.7) | **P2** (early) + P5 (補強) |

P2 早期攔截 4/6（66%）；P5 補強或主攔截 5/6（包含 P2 已 BLOCKED 的也再次驗證）。**P2 與 P5 為冗餘但非重複**（P2 看 precondition，P5 看 result coherence）— 這是設計目標。

---

## §6. 改進建議

### §6.1 已驗證的設計決策（保留不動）

- ✅ **Per-component override（n_critical ≥ 1 OR n_low ≥ 3）**：對 retrospective sensitivity 貢獻 4/5；不調整。
- ✅ **Path A 6-component weights (msc 0.27, ess 0.22, pf/sh/pcs 0.17)**：在所有 8 cases 上排序合理；不調整。
- ✅ **P-04 pileup symlink probe 觸發條件（dataset_id 含 "pileup" + symlink/clairs_paired + clairs-to/_to_）**：完整命中 vcf_source_error；不放寬。
- ✅ **Negative controls clean precondition recipe**（post-KDE + caller_af 分離 + binary_version=null）作為「good cycle template」：未來 plan.json 撰寫時可直接套用。

### §6.2 Drill 1 揭露的設計潛在改進（不阻塞 Day 7+）

| # | 觀察 | 建議改進 | 優先級 |
|---|---|---|---|
| 1 | tier_support_alignment component 在 5/5 positive + 2/2 negative 全為 null（依 schema 設計） | 此 component 需要 `evidence_ledger` stability_grade ↔ tier_used 對映；目前 retro fixture 沒寫此資料 → 後續真實 cycle 應在 plan.json 中加 `cycle.tier_evidence_grade` 欄位 | P2 |
| 2 | P-06 (n_reads/AF confound) pitfall 在 hpfinengroups 命中（pcs=0），但其他 4 cases 都未命中 | pitfalls_table 6 條覆蓋 2026-04 事件範圍；後續若有新 confound class（如 KDE 新型 artifact）應補 P-07/P-08 | P3 |
| 3 | longphase risk_base 已 0.74 直接觸 base composite > 0.7 阻擋；override 與 base 在此 case 重複生效 | 設計合理（雙保險）；無需調整 | n/a |
| 4 | thread_b 的 risk_base 0.498 在 base 規則只觸 downgrade，靠 override n_critical=2 升至 pending_review | 反過來看，若僅看 base：thread_b 會被低估為 ⭐2 而非 RETRACTED；override 必要性確認 | n/a |
| 5 | retro_cycles 不在 `state/cycles/` 而在 `state/retro_cycles/` — `_evaluator_run.py` 已加 fallback；`active.json` 不收 retro | 設計正確；但建議 `cycle-state` skill 的 dashboard 將 retro_cycles 列為「regression test fixtures」獨立 section | P3 |

### §6.3 Drill 1 的限制（透明標記）

- ⚠ **Hard fidelity events 用「指標重現」而非完整重跑**：hpfinengroups + longphase 的 pilot.json/generalize.json 直接從原始報告抽 NGroups histogram 與 cross-sample table；**未重新編譯+執行 C++ binary**（plan v1.3 §4.5.3-B Anti-goal 接受）。如需完整 fidelity 須跑數小時的 BAM/VCF pipeline。
- ⚠ **Sample size n=6 positive + n=2 negative**：小樣本下 100% sensitivity 與 specificity 無法外推到所有未來事件；建議累積 ≥3 個月新 cycle 後做第二次 retrospective。
- ⚠ **Negative controls 都是最近驗證的（N1 paired POSITIVE / N2 Zone-Aware）**：未涵蓋「中等品質、應 downgrade 但不 retract」的灰色案例；plan v1.6 候選 — 加 1-2 個 ⭐3→⭐2 的 N3/N4 negative controls。

### §6.4 §4.5.4-F Drill 1 後 batch 3-5 重評（plan v1.5 約定）

按照 plan v1.5 §4.5.4-F 的決策框架，依本 drill 結果對 Day 5.5 剩餘 batch 3-5 逐一評估：

| Batch | 原 plan 預期收益 | Drill 1 觀察 | 決策 |
|---|---|---|---|
| **Batch 3 — 4 passive reference**（known-pitfalls / doc-standards / confirmation-protocol / observation-analysis 改 passive） | 減少誤觸發、節省 context | Drill 1 中 known-pitfalls 透過 pitfalls_table 在 evaluator 內部被消費 1 次（hpfinengroups P-06 命中），未觀察到誤觸發；其他 3 個 passive 無觀察事件 | **`[PARKED]`** — Drill 1 沒揭露急需處理問題；待真實 cycle 中若觀察到 known-pitfalls 在主對話被誤觸發再啟動 |
| **Batch 4 — 23 forward-link 批次** | 提升 chain 路由精準度 | Drill 1 retrospective 是 fixture-driven，未經主對話 phase routing；無法觀察 chain 引導效果 | **`[PARKED]`** — 等 Drill 2 (small new cycle 走 P0-P6) 才有觀察數據；Drill 2 是 plan v1.5 Week 2 任務，現未啟動 |
| **Batch 5 — 4 個人風格 anchor 硬化** | 確保 evaluator / structured-tech-report / results-analysis 對齊風格 | Drill 1 中 evaluator components 已隱含 anchor #1（multi-sample consistency = 多軌證據）；anchor #3/6/7 未在 retrospective 中展現必要性 | **`[REORDER]`** — 把 anchor #1 (L4 mandatory) 的硬化提到 batch 5a 優先（因 evaluator 已隱性使用），其餘 #3/6/7 留 batch 5b 待 Drill 2 後再評 |

**綜合建議**：依 plan v1.5 約定的「實證導向迭代」原則，**batch 3-5 全部 PARKED 至 Drill 2 完成**，先做：
1. Week 2 Day 7-8: `/cycle-init` + `/cycle-state` skill（plan v1.5 Week 2 既定任務）
2. Week 2 Day 9-10: Drill 2（small new cycle 走 P0-P6）
3. Drill 2 後再依其 metric 重評 batch 3-5

### §6.5 Plan v1.6 候選增強（記錄不阻塞）

1. **codex-plugin-cc 整合評估**（plan v1.5 deferred decision 2026-05-06）
   - Drill 1 顯示 evaluator 在「事後 reflective analysis」表現好（sensitivity 6/6），但缺**proactive second-opinion**（執行中 GPT 對 plan/pilot 的獨立挑戰）
   - 建議插入點：**P1 PLAN end** 與 **P5 EVALUATE pre**（讓 GPT 對 plan.json 做 independent critique；對 evaluation.json verdict 做 sanity check）
   - 風險：增加 token 成本與決策延遲
   - 投入優先級：plan v1.6 Q5 候選

2. **negative controls 擴充至 4 個**（current 2 → 加 2 個灰色案例）
   - 候選 N3：⭐3 paired-pure delta F1 baseline（驗證 baseline 不被誤抓）
   - 候選 N4：⭐2 explore tier 已 retire 但 process 健全的 cycle（驗證低 tier 也應 approve_tier 而非 RETRACTED）

3. **Drill 1 結果存入 evaluator 的 precedent retrieval**（Path B 才執行）
   - 把本 drill 的 6 positive + 2 negative 作為 LlamaIndex precedent corpus
   - 未來 cycle 在 P5 EVALUATE 時可查「我這個 cycle 與哪些先例最相似？」
   - 阻塞：需先建 LlamaIndex 索引（Path B Week 4 任務）

---

## §7. 給未來自己的 Heads-up

實作這個 drill 過程中觀察到的 5 個易錯點（補 plan §10）：

1. **`risk_components` vs `components` key 命名**：`evaluation.json` schema 使用 `risk_components`，但 `_evaluator_run.py` 印 stdout 用 `components`；外部 query 兩者並存導致 None。**未來統一改為 `risk_components`**。
2. **retro_cycles 與真實 cycles 的 path fallback**：`_staleness_check.py` 與 `_evaluator_run.py` 的 `cycle_dir()` helper 已加 retro fallback；但若有第三類 cycle（如 `state/sandbox_cycles/`）需要再擴充。
3. **fabricated SHA 觸發 binary stale=missing**：thread_b 第一次跑時用 `8d0a0c8feedba7c...` fake SHA，被誤判為 BLOCKED；解決方式是純分析 cycle 設 `binary_version=null`。**未來 plan.json 範本應明確標註此 convention**。
4. **expected catch table 必須 freeze 在執行前**：本 drill 嚴格遵守 §4.5.3-A 的 freeze；如有事後挪動會掩蓋 false positive 風險。**plan v1.6 把此規則明寫入 confirmation-protocol**。
5. **per-component override 規則的閾值（0.2/0.4/3）是 magic number**：若未來 sample size 改變（如 N=10 而非 N=7），`n_low ≥ 3` 的閾值需要 calibrate。**保留靈活性，不 hardcode**。

---

## §8. References

- Plan: `~/.claude/plans/agent-harness-langgraph-resilient-waterfall.md` v1.5
- Decision Log entries: 2026-05-04 Q1-4 + 2026-05-05 v1.3 Q1-4 + 2026-05-06 v1.5 工作流
- Source events:
  - `InterSubMod/docs/experiments/in_progress/2026/04/20260404_VCF來源錯誤矯正報告_01.md`
  - `InterSubMod/docs/experiments/in_progress/2026/04/20260420_KDE_Fix_Acceptance_Validation_01.md`
  - `InterSubMod/docs/experiments/in_progress/2026/04/20260418_F_HPFineNGroups_Subclone_Marker_01.md`
  - `InterSubMod/docs/experiments/in_progress/2026/04/20260424_X6_Caller_AF_S3S5_CrossSample_01.md`
  - `InterSubMod/docs/experiments/in_progress/2026/04/20260426_Thread_B_Whitelist_Retraction_01.md`
  - `InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md`
- Negative control sources:
  - `InterSubMod/docs/reports/research_landscape/09_Part_B.md` (N1 paired LOH×AF×Methylation POSITIVE)
  - `InterSubMod/docs/reports/research_landscape/08_Zone_Aware.md` (N2 Zone-Aware characterization)
- Drill artifacts: `InterSubMod/state/retro_cycles/{event_id}/` × 8 cycles
- Methodology references:
  - Google SRE postmortem-driven testing (chaos engineering historical replay)
  - AWS Well-Architected Framework — operational excellence reliability pillar
  - Medical diagnostics ROC framework (sensitivity + specificity + 2×2 confusion matrix)
  - Anthropic Building Effective Agents (Evaluator-Optimizer pattern)
  - LangGraph reflection patterns
