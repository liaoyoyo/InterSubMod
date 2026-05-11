---
id: ism-kb-09-conclusions-negative-20-paired-f1-filter
name: "N20: Paired F1-filter 方向放棄（2026-04-21 決策）"
description: "2026-04-21 HCC1395 R1 pilot paired arm Fisher_Frac_Sig raw AUC 0.726、殘差化 0.698 [CI 0.534, 0.831] CI 跨隨機；F pilot subset FP=10 不足；cross-mode TP 99.5% 飽和；使用者決定放棄 paired F1-filter，保留 characterization-only（須全域 region）。"
status: archived
last_verified: 2026-04-21
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_paired_f1_filter_abandoned.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, paired, fisher-frac-sig, f1-filter, characterization]
canonical_paths: [09_conclusions/negative/N20_paired_f1_filter.md]
alias_paths: []
---

# N20: Paired F1-filter 方向放棄

## Status
`concluded-negative` for F1-filter；characterization-only 保留（須跳出 F pilot subset 用全域 region 才有意義）。

## Context
HCC1395 R1 pilot paired arm 的 `Fisher_Frac_Sig` raw AUC 0.726 看似突破 CL-008 Beyond-AUC ≤0.58 ceiling，疑似為 paired 模式的 F1-filter 候選。需殘差化與跨樣本驗證。

## Decision-tested
Paired `Fisher_Frac_Sig` 是否為穩定獨立的 F1-filter 候選（殘差化 AUC CI 不跨隨機 + 跨樣本一致）。

## Method
- Paired arm Fisher_Frac_Sig raw + residualized AUC
- Bootstrap CI
- Cross-mode concordance（F pilot subset 內 paired vs TO 雙陽性率）
- Paired Fisher_Frac_Sig vs TO SampleASM_Delta 相關性

## Result
- **Raw AUC**：0.726（疑似突破 ≤0.58 ceiling）
- **殘差化 AUC**：0.698 [CI **0.534**, 0.831] — **CI 下界跨入隨機**
- F pilot subset FP=10 過小（bootstrap 不穩）
- **Cross-mode concordance**：F pilot subset 內 TP 99.5% 飽和（both_TP 439/441），paired arm 在 subset 內已無 F1 空間
- **獨立性**：paired Fisher_Frac_Sig 與 TO SampleASM_Delta ρ=−0.162（獨立）→ characterization 有獨立基礎但 F1-filter 無效

## Consequences
- **根因**：F pilot subset 在 paired arm 已達 TP 99.5% 飽和（no room for F1-filter）；殘差化後 CI 跨隨機說明原 raw AUC 為 confound-driven。
- **替代方向**：
  - Paired Fisher_Frac_Sig 不再視為 F1-filter 候選，除非（1）R12 跨樣本在某樣本發現大量 FP region 且（2）paired Fisher_Frac_Sig 在該樣本獨立於 TO SampleASM_Delta 且 AUC ≥ 0.65
  - Characterization-only 使用時**必須跳出 F pilot subset（全域 region）**
  - TO arm `SampleASM_Delta` 仍是 P0 F1-filter 主候選（CL-025a ⭐3，殘差化 0.610 穩定）— 但後已被 N02 R1-Global 否決
- **決策**：2026-04-21 使用者明確指示「放棄 paired F1-filter 方向，轉 characterization-only」。Registry CL-025b status 更新為 "concluded (abandoned for F1-filter, characterization-only per user 2026-04-21)"。

## References
- MEMORY：`project_paired_f1_filter_abandoned.md`
- `research/F_hpfinengroups_deepening/observations/step7_crossmode_concordance.md`
- `research/F_hpfinengroups_deepening/observations/step7_hcc1395_normal_paired_pilot.md`
- Registry：`docs/reports/research_landscape/10_Research_Chain_Registry.md` §CL-025b / CL-025c
