---
postmortem_id: 20260518_readparser_germline_hp_only_phase1
topic: ReadParser --germline-hp-only flag Phase 1
verdict: CONDITIONAL NEGATIVE (filter-level) → REOPENED via V6 (characterization-level POSITIVE)
date: 2026-05-18
framing: blameless (machinery, not individual)
related_cycles:
  - phase1_readparser_germline_hp_only (2026-04-21, commit 775027d)
  - v6_priority_bug_audit (2026-05-10 to 2026-05-12)
next_recall: 2026-06-17
next_recall_2: 2026-08-17
---

# Postmortem: ReadParser `--germline-hp-only` Phase 1 NEGATIVE → V6 Reopen

> **使用 `/scientific-rigor` §9.2 inline template**；同步 §8.3 Reflexion buffer + §8.3.1 Reopen threshold 評估。

## §1 Summary（≤3 行）

- ReadParser `--germline-hp-only` flag 機制完全正確（chr19 smoke + 12 unit tests pass），但 0421 Phase 1 HCC1395 TO 全量 **無 TSV 特徵 AUC delta ≥ +0.02** → filter-level NEGATIVE。
- 本 NEGATIVE 觸發 HPFineNGroups marker 機制重新詮釋（從 methylation subclone marker → LOH-constrained phasing signature），引出新主軸候選。
- 0510-0511 V6 binary patch 在 chr19 + 4 樣本 cross-sample 證實 marker 真實生物學訊號（schema collapse 後物理屬性仍可區辨 TP/FP）→ **characterization-level reopen 成立**；filter-level NEGATIVE 維持不撤回。

## §2 Timeline

| Date | Event | Artifact |
|------|-------|---------|
| 2026-04-21 | Phase 0 chr19 smoke（615 TP）+ 12 unit tests pass；機制獨立、demotion 數學守恆 | commit `775027d` (`refactor/phase1-safety`) |
| 2026-04-21 | Phase 1 HCC1395 TO 全量（28,509 TP + 11,606 FP）；HPFineNGroups Δ=-0.026、HPMergedDelta -0.025、NHP3 -0.035 | `InterSubMod/docs/experiments/in_progress/2026/04/20260421_ReadParser_GermlineHPOnly_HCC1395_01.md` |
| 2026-04-21 | flag=on 時 NG≥3 完全消失（0/28,509 TP regions）→ 原 ⭐4 subclone marker 結論觸發 audit | 同上 |
| 2026-04-23 | 週報定案：HPFineNGroups 機制從 methylation subclone marker → LOH-constrained phasing signature；marker 降級 ⭐4→⭐3 pipeline-dependent | `project_hpfinengroups_subclone_marker.md` |
| 2026-05-10 | V6 binary patch（移除 `HaplotagProcess.cpp:537-548` Layer 1.5 else if 分支）+ chr19 cross-tab 768 regions：marker rate 94.7% + NG_on=2 cell 91.5% → chr19 conditional POSITIVE | `InterSubMod/research/paired_priority_bug_audit/05_V6C_phaseB_findings.md` |
| 2026-05-11 | V6 Phase D 4 樣本 cross-sample 一致 POSITIVE；hp=1-1:hp=2-1 ratio 平均 0.93（vs HCC1395 V5 1.86, baseline 17.3） | `InterSubMod/research/paired_priority_bug_audit/08_phaseD_v6_cross_sample_findings.md` |
| 2026-05-12 | V5 vs V6 SP1/2/3 反向 trade-off 量化：V5 在 phasing-weak 區、V6 在 germline-absent 區各勝 | `project_v5_v6_tradeoff_sp123.md` |

## §3 Root Cause（machinery, not blame）

**多層次根因，非單一錯誤**：

### 3.1 執行層
Filter-level evaluation（AUC delta on TP/FP discrimination）對「機制層 schema collapse 的數學必然性」不敏感。flag=on 後 NG≥3 = 0 是 bucket demote（HP:i:11/21/33 → 0）的數學必然，不能反推為 marker artifact 證據 — 但 AUC delta 看不出這個區別。

### 3.2 方法論層
原 Plan baseline 0.836（引自 `02_Self_Phasing根因.md`）vs 實測 V3-Fixed TO BAM = 0.549，差距 0.287 未在 Phase 0 後先 reconcile。投入 Phase 1 全量前缺 §7.2 reproducibility checkpoint（dataset 版本 / haplotag 版本 / VCF 子集 mismatch 來源不明）。

### 3.3 概念層
HPFineNGroups 命名 + 既有文獻引用讓 reader 預設「subclone methylation marker」，掩蓋了它本質是 4-bucket HP × variant occupancy count（非 methylation bimodality）。當 flag toggle 改變 bucket schema，marker 必然失效，但這是 **resolution loss**，不是 **signal disappearance**。

## §4 What Went Well

1. **Phase 0 chr19 smoke 先跑**（615 TP），未直接全量 → 機制正確性早期建立；12 unit tests + audit 獨立保住 demotion 數學守恆。
2. **願意質疑既有 ⭐4 結論不護航**：發現 flag=on NG≥3=0 立即列入「marker 可能撤回」候選，未因投入沉沒成本維護原結論。
3. **NEGATIVE 觸發深度根因調查 → 新主軸誕生**：LOH-constrained phasing discovery（4-22）+ V6 binary patch（5-10）皆是本 NEGATIVE 的下游 productive failure 產物。
4. **V6 patch 完整三向驗證**：chr19（局部）+ 全基因組（intra-sample）+ 4 樣本（cross-sample）符合 `/scientific-rigor §5` 多方驗證原則。
5. **修正保留為 opt-in（default=off）**：未因 NEGATIVE 直接 revert，給研究者保留乾淨 HP 分群選項。

## §5 What Went Poorly

1. **Plan baseline reconcile 缺前置 staleness check**：應在 P2 PRECHECK 階段攔下 0.836 vs 0.549 落差，避免 Phase 1 投入全量後才事後標註。
2. **HPFineNGroups 原 ⭐4 升級時未要求 flag-toggle 重驗證**：marker 機制依賴單一 binary 行為（V3F 含 priority bug），升級驗證未覆蓋 binary parameter sweep。
3. **未先做 TP/FP HP-tag density disparity check**：TP 27.4 tags/site vs FP 27.7 tags/site 幾乎相同 — 若先驗證此前置，可預測 filter 只能去 noise 不能產生 signal，避免投入全量 Phase 1 AUC 評估。
4. **`/known-pitfalls` 缺對應陷阱條目**：Pipeline-dependent marker 重驗證原則未文件化，導致重複可能性。

## §6 Action Items

| # | 項 | Owner | Due | Status |
|---|----|------|-----|-------|
| A1 | `/known-pitfalls` 新增 P-XX「Pipeline-dependent marker 必 flag-toggle 重驗證」 | AI | 下輪 docs round | open |
| A2 | V6 Phase C 7 樣本 V5 BAM 全量 ISM × 2 flag 擴展（joint NG ∧ AF<0.4 ∧ NR≥80 ∧ NonLOH 完整 filter） | researcher | TBD | open |
| A3 | HPFineNGroups marker 機制描述全域 search & replace（subclone → phasing signature） | AI | 下輪 docs round | open |
| A4 | Plan baseline reconcile 流程加入 `/check-staleness` P2 PRECHECK template | AI | infra round | open |
| A5 | TP/FP feature-density disparity pre-check 加入 `/methodology-audit` Step 1 量化現況 checklist | AI | 下輪 methodology round | open |

## §7 Reopen Threshold 評估（依 `/scientific-rigor` §8.3.1）

| 條件 | 是否滿足 | 證據 |
|------|--------|------|
| **C1 新數據** | ⚠ 部分 | V6 haplotag + V5 phasing 4 樣本 cross-sample data；尚缺 7 樣本完整擴展 |
| **C2 新方法** | ✅ 滿足 | V6 binary patch 是 systematic 新框架（移除 Layer 1.5 else if 繞過 priority bug confound） |
| **C3 新前置** | ✅ 滿足 | Priority bug audit（hp=1-1:hp=2-1 ratio 17.3 → 0.93）解除原 confound 機制 |

**結論**：
- **Filter-level NEGATIVE 維持不撤回**（AUC delta 在 V6 下仍未測，但 §5 What Went Poorly #3 預測仍不會翻盤）
- **Characterization-level reopen 已落地**（0510-0511 V6 Phase B/D conditional POSITIVE）
- **新 cycle 標籤**：`reopen:phase1_readparser_germline_hp_only` → 待 V6 Phase C 7 樣本完整擴展

## §8 Reflection Buffer（給下次同方向研究的 agent）

```markdown
## REFLECTION

**警示指標**：未來提案「filter 機制修正提升 AUC」前，先跑 mid-AUC feature inventory；若該 feature group 內無 mid-AUC（>0.55）HP-related feature，AUC delta 不會翻盤 — 應改 characterization-level 評估而非 filter-level。

**根因**：somatic tag 在 TP/FP 上密度幾乎相同（TP 27.4/site vs FP 27.7/site），filter 只能去 noise 不能產生 signal；marker 機制依賴 bucket schema，flag toggle 改變 schema 後 marker 必失效但不代表 signal 消失。

**改進方向**（若要重試）：
1. 必先做 TP/FP HP-tag density disparity check（prerequisite）
2. Marker 在 binary parameter sweep 下重驗證（≥2 個 flag 組合）
3. 評估「characterization-level (region property)」而非「filter-level (AUC delta)」
4. Plan baseline vs 實測落差先 reconcile（避免 sunk cost）

**Spaced recall dates**：2026-06-17（30d）+ 2026-08-17（90d）

**Productive failure 價值**：此 NEGATIVE 觸發 LOH-constrained phasing discovery + V6 binary patch + priority bug audit 三條下游主軸；高轉化價值，非浪費。
```

## §9 Links

### Reports
- 原 NEGATIVE Phase 1 report: `InterSubMod/docs/experiments/in_progress/2026/04/20260421_ReadParser_GermlineHPOnly_HCC1395_01.md`
- Smoke report: `InterSubMod/docs/experiments/in_progress/2026/04/20260421_ReadParser_GermlineHPOnly_smoke_01.md`
- V6 binary validation: `InterSubMod/research/paired_priority_bug_audit/07_V6_validation_findings.md`
- V6 Phase B chr19: `InterSubMod/research/paired_priority_bug_audit/05_V6C_phaseB_findings.md`
- V6 Phase D 4 樣本: `InterSubMod/research/paired_priority_bug_audit/08_phaseD_v6_cross_sample_findings.md`
- V6 production documentation: `InterSubMod/docs/reports/validated/2026/05/20260511_V6_binary_complete_documentation_01.md`
- HP tag priority bug engineering report: `InterSubMod/docs/reports/validated/2026/05/20260514_HP_tag_priority_bug_engineering_report_01.md`

### Memory entries
- `project_readparser_germline_hp_only_phase1_negative.md` — 原 NEGATIVE 結論
- `project_hpfinengroups_subclone_marker.md` — ⭐4→⭐3 降級 + V6 reopen
- `project_loh_subclone_af_methylation_positive.md` — TO 機制 pivot
- `project_loh_constrained_phasing_discovery.md` — 新主軸候選
- `project_v5_v6_tradeoff_sp123.md` — V5 vs V6 trade-off

### Commits
- Phase 0 + smoke: `775027d` on `refactor/phase1-safety`
- V6 binary patch: 詳見 priority bug audit 系列 commit

### Cross-references
- `/scientific-rigor` §8.3.1 Productive Failure + Reopen Threshold（C2 + C3 觸發本案 reopen）
- `/known-pitfalls`（A1 待新增 P-XX「Pipeline-dependent marker」）
- `InterSubMod/.claude/CLAUDE.md` §1 Hard Gate（NO-GO 判定不可事後改寫 — 本案 filter-level NO-GO 維持不撤回，遵守此規則）

---

**Status**: postmortem complete; next spaced recall 2026-06-17；待 V6 Phase C 7 樣本擴展 cycle 後本檔加 §10 update。
