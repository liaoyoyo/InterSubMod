---
name: research-context-loader
description: Loads InterSubMod research landscape (10 reports + experiments INDEX + concepts) on-demand when user starts research analysis, hypothesis verification, or asks about prior conclusions. USE WHEN 「延續研究」「研究方向」「結論驗證」「landscape」「過去結論」「研究歷史」「inj research」「pivot」「P0-P5 cycle」「假說 cycle」. SKIP WHEN simple build/test, code-only edits, doc writing, single-file analysis, no research context needed.
type: skill
allowed-tools: Read, Grep, Glob
---

# Research Context Loader

## Phase & Chain Position

Cross-cutting context-provider. Activates BEFORE research analysis to load landscape context, separate from the live focus snapshot in CURRENT_FOCUS.md (which is always-loaded via CLAUDE.md).

## Dependencies

- **Reads**: 10 landscape reports + experiments INDEX + concepts (paths in playbook.md)
- **Used by**: Any research analysis (research-loop, feature-layered-observation, conclude-research, pivot-direction)
- **Replaces**: Old CLAUDE.md "6 必讀清單" landscape rows (now conditional)

## Failure Mode & Diagnostics

| Symptom | Cause | Diagnosis |
|---------|-------|-----------|
| User asks about a past conclusion, model gives wrong/outdated info | Landscape not loaded for this session | Re-invoke this skill before answering |
| Token usage spikes after research session start | Loaded entire landscape unnecessarily | Use Tier 2 (just INDEX) instead of Tier 3 (full files) |
| Conflicting conclusions in same session | Loaded subset of landscape, missed key file | Load `00_INDEX` + `06_結論穩定性審查` first to discover relevant deeper files |

## 3-Tier Loading Strategy

See `playbook.md` for the decision tree. Brief version:

- **Tier 1** (always-loaded by CLAUDE.md): `docs/CURRENT_FOCUS.md` only
- **Tier 2** (this skill, light): `docs/experiments/INDEX.md` + `docs/reports/research_landscape/00_INDEX.md` + `docs/concepts/2026/04/20260409_研究構想總索引_01.md` (~3-5k tokens, common research question)
- **Tier 3** (this skill, deep): specific landscape files based on question topic (~20-25k tokens, deep dive question)

## Quick reference (from CLAUDE.md table)

| 我想知道... | 去哪裡找 |
|-------------|---------|
| TO FP 問題全貌 | `docs/reports/research_landscape/01_TO_FP問題全貌.md` |
| Self-phasing 根因與影響 | `docs/reports/research_landscape/02_Self_Phasing根因.md` |
| 哪些 ISM 特徵可信 | `docs/reports/research_landscape/03_ISM分析價值界定.md` |
| 暫停判定與修正後預期 | `docs/reports/research_landscape/04_暫停判定與重評估.md` |
| 完整證據鏈推論 | `docs/reports/research_landscape/05_證據鏈總覽.md` |
| 結論穩定性評分 | `docs/reports/research_landscape/06_結論穩定性審查.md` |
| LOH/CN/AF 三維度統合 | `docs/reports/research_landscape/07_LOH_CN_AF_研究總整理.md` |
| Zone-Aware Framework 驗證歷程 | `docs/reports/research_landscape/08_Zone_Aware.md` |
| Part B 質疑驗證（HPFineNGroups 升級） | `docs/reports/research_landscape/09_Part_B.md` |
| 研究鏈完整登記簿 | `docs/reports/research_landscape/10_Research_Chain_Registry.md` |
| Truth set 答案與 F1 比對協議（7 樣本） | `docs/data_specs/20260422_truth_sets_and_f1_protocol_01.md` |
| 啟動壓縮上下文、任務順序、待確認問題 | `docs/references/manual/20260424_AI啟動壓縮上下文與研究索引_01.md` |
