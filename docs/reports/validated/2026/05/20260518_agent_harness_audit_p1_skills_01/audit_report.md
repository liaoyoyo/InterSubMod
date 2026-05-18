# P1 Skills Audit — InterSubMod Agent Harness

**Audit date**: 2026-05-18 14:30
**Skills scanned**: 42
**Avg line count**: 188 / **Avg desc len**: 432 chars

## Executive Summary

- 🔴 **CRITICAL — YAML invalid**: 5 skills（description 以 `**` 開頭被 YAML 解讀為 alias reference）
- 🟠 **D2 內容診斷段缺**: 22 / 42 skills 缺 Phase chain / Dependencies / Failure Mode 段
- 🟡 **D3 cross-ref 缺**: 17 skills 無與 /scientific-rigor 連結（含 5 個非工具類真問題）
- ✅ **D1/D4/D5 全綠**: frontmatter / trigger words / staleness 都 OK

## 6 維度評分

| Dim | ✅ | ⚠ | ❌ | 維度說明 |
|-----|----|----|----|---------|
| **D1** | 42 | 0 | 0 | frontmatter 完整度 |
| **D2** | 14 | 6 | 22 | 內容診斷段 |
| **D3** | 11 | 14 | 17 | cross-ref to /scientific-rigor |
| **D4** | 42 | 0 | 0 | trigger 詞 USE WHEN/SKIP WHEN |
| **D5** | 42 | 0 | 0 | staleness mtime |

## 詳細表 — 42 skills × 6 維度

| Skill | Group | D1 frontmatter | D2 診斷 | D3 cross-ref | D4 trigger | D5 mtime | YAML | Lines | Age (d) |
|-------|-------|----------------|---------|--------------|------------|----------|------|-------|---------|
| `confirmation-protocol` | 元方法論 | ✅ | ❌ | ✅ | ✅ | ✅ | ✅ | 288 | 0 |
| `cycle-state` | 元方法論 | ✅ | ✅ | ⚠ | ✅ | ✅ | ✅ | 164 | 7 |
| `fast-learning-coach` | 元方法論 | ✅ | ❌ | ✅ | ✅ | ✅ | ✅ | 281 | 0 |
| `known-pitfalls` | 元方法論 | ✅ | ❌ | ✅ | ✅ | ✅ | ✅ | 180 | 0 |
| `research-context-loader` | 元方法論 | ✅ | ✅ | ❌ | ✅ | ✅ | ✅ | 51 | 7 |
| `scientific-rigor` | 元方法論 | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | 599 | 0 |
| `check-staleness` | 7-Phase Waterfall | ✅ | ❌ | ✅ | ✅ | ✅ | ✅ | 115 | 0 |
| `conclude-research` | 7-Phase Waterfall | ✅ | ❌ | ⚠ | ✅ | ✅ | ✅ | 231 | 7 |
| `cycle-init` | 7-Phase Waterfall | ✅ | ✅ | ⚠ | ✅ | ✅ | ✅ | 174 | 7 |
| `feature-layered-observation` | 7-Phase Waterfall | ✅ | ✅ | ⚠ | ✅ | ✅ | 🔴 | 326 | 7 |
| `multi-sample-consistency` | 7-Phase Waterfall | ✅ | ❌ | ⚠ | ✅ | ✅ | ✅ | 198 | 7 |
| `research-loop` | 7-Phase Waterfall | ✅ | ✅ | ⚠ | ✅ | ✅ | 🔴 | 442 | 7 |
| `run-evaluator` | 7-Phase Waterfall | ✅ | ⚠ | ⚠ | ✅ | ✅ | ✅ | 184 | 7 |
| `cpp-change` | 程式開發 | ✅ | ⚠ | ❌ | ✅ | ✅ | ✅ | 267 | 0 |
| `infra-ops` | 程式開發 | ✅ | ✅ | ❌ | ✅ | ✅ | ✅ | 55 | 7 |
| `methodology-audit` | 程式開發 | ✅ | ❌ | ✅ | ✅ | ✅ | ✅ | 136 | 0 |
| `verification-loop` | 程式開發 | ✅ | ✅ | ✅ | ✅ | ✅ | 🔴 | 212 | 0 |
| `citation-verification` | 文件管理 | ✅ | ⚠ | ❌ | ✅ | ✅ | ✅ | 48 | 7 |
| `data-audit` | 文件管理 | ✅ | ❌ | ❌ | ✅ | ✅ | ✅ | 162 | 7 |
| `doc-standards` | 文件管理 | ✅ | ❌ | ❌ | ✅ | ✅ | ✅ | 107 | 9 |
| `memory-consolidation` | 文件管理 | ✅ | ❌ | ✅ | ✅ | ✅ | ✅ | 135 | 0 |
| `html-preview` | 報告生成 | ✅ | ✅ | ❌ | ✅ | ✅ | ✅ | 137 | 4 |
| `html-report-build` | 報告生成 | ✅ | ✅ | ❌ | ✅ | ✅ | ✅ | 269 | 4 |
| `myPPT` | 報告生成 | ✅ | ❌ | ⚠ | ✅ | ✅ | ✅ | 75 | 7 |
| `pptx-build` | 報告生成 | ✅ | ❌ | ⚠ | ✅ | ✅ | ✅ | 229 | 7 |
| `report` | 報告生成 | ✅ | ❌ | ❌ | ✅ | ✅ | ✅ | 99 | 7 |
| `results-report` | 報告生成 | ✅ | ❌ | ❌ | ✅ | ✅ | ✅ | 181 | 7 |
| `structured-tech-report` | 報告生成 | ✅ | ❌ | ⚠ | ✅ | ✅ | ✅ | 201 | 7 |
| `weekly-report` | 報告生成 | ✅ | ❌ | ❌ | ✅ | ✅ | ✅ | 246 | 9 |
| `auc-confound-guard` | 研究專用 | ✅ | ❌ | ✅ | ✅ | ✅ | ✅ | 162 | 0 |
| `image-gen` | 研究專用 | ✅ | ✅ | ❌ | ✅ | ✅ | ✅ | 82 | 7 |
| `image-vision-check` | 研究專用 | ✅ | ✅ | ❌ | ✅ | ✅ | ✅ | 72 | 7 |
| `init-research` | 研究專用 | ✅ | ⚠ | ⚠ | ✅ | ✅ | ✅ | 144 | 7 |
| `inject-hypothesis` | 研究專用 | ✅ | ❌ | ❌ | ✅ | ✅ | ✅ | 150 | 7 |
| `observation-analysis` | 研究專用 | ✅ | ❌ | ❌ | ✅ | ✅ | ✅ | 139 | 7 |
| `pivot-direction` | 研究專用 | ✅ | ⚠ | ⚠ | ✅ | ✅ | ✅ | 99 | 7 |
| `problem-framing-ideation` | 研究專用 | ✅ | ✅ | ⚠ | ✅ | ✅ | 🔴 | 269 | 7 |
| `provenance-tier-audit` | 研究專用 | ✅ | ✅ | ✅ | ✅ | ✅ | 🔴 | 224 | 0 |
| `research-dashboard` | 研究專用 | ✅ | ❌ | ⚠ | ✅ | ✅ | ✅ | 98 | 7 |
| `results-analysis` | 研究專用 | ✅ | ❌ | ❌ | ✅ | ✅ | ✅ | 210 | 7 |
| `review-evidence` | 研究專用 | ✅ | ⚠ | ❌ | ✅ | ✅ | ✅ | 200 | 7 |
| `validation-protocol` | 研究專用 | ✅ | ❌ | ✅ | ✅ | ✅ | ✅ | 270 | 0 |

## 🔴 Critical Fix List (5 YAML invalid)

這 5 個 SKILL.md description 以 `**` 開頭，違反 YAML 規範（`**` 視為 alias reference）。Anthropic skill loader 寬鬆度未知，是 latent bug。

| Skill | 修法 |
|-------|------|
| `feature-layered-observation` | description 用 `"..."` quoted string 包起來，或移除開頭 `**` markdown |
| `problem-framing-ideation` | description 用 `"..."` quoted string 包起來，或移除開頭 `**` markdown |
| `provenance-tier-audit` | description 用 `"..."` quoted string 包起來，或移除開頭 `**` markdown |
| `research-loop` | description 用 `"..."` quoted string 包起來，或移除開頭 `**` markdown |
| `verification-loop` | description 用 `"..."` quoted string 包起來，或移除開頭 `**` markdown |

**修法範例**：
```yaml
# ❌ Before
description: **P3 PILOT 主要分析 skill** — ...
# ✅ After (quoted)
description: "**P3 PILOT 主要分析 skill** — ..."
# ✅ After (no markdown)
description: P3 PILOT 主要分析 skill — ...
```

## 🟠 D2 內容診斷段缺 (Priority by group)

依 user feedback `feedback_skill_md_must_state_dependencies_and_diagnostics`，所有 SKILL.md 應含 3 段：
- `Phase & Chain Position`（在工作流哪一階段）
- `Dependencies` (Uses / Used-by / Reads / Writes)
- `Failure Mode & Diagnostics`（失敗訊號 + 排查路徑）

**非工具類 13 skill 缺診斷段（高優先）**：

| Skill | Group | Lines | 缺項 |
|-------|-------|-------|------|
| `confirmation-protocol` | 元方法論 | 288 | Phase chain, Deps, Failure |
| `fast-learning-coach` | 元方法論 | 281 | Phase chain, Deps, Failure |
| `known-pitfalls` | 元方法論 | 180 | Phase chain, Deps, Failure |
| `check-staleness` | 7-Phase Waterfall | 115 | Phase chain, Deps, Failure |
| `conclude-research` | 7-Phase Waterfall | 231 | Phase chain, Deps, Failure |
| `multi-sample-consistency` | 7-Phase Waterfall | 198 | Phase chain, Deps, Failure |
| `methodology-audit` | 程式開發 | 136 | Phase chain, Deps, Failure |
| `memory-consolidation` | 文件管理 | 135 | Phase chain, Deps, Failure |
| `auc-confound-guard` | 研究專用 | 162 | Phase chain, Deps, Failure |
| `inject-hypothesis` | 研究專用 | 150 | Phase chain, Deps, Failure |
| `observation-analysis` | 研究專用 | 139 | Phase chain, Deps, Failure |
| `results-analysis` | 研究專用 | 210 | Phase chain, Deps, Failure |
| `validation-protocol` | 研究專用 | 270 | Phase chain, Deps, Failure |

## 🟡 D3 cross-ref 缺 — 5 個應補（非工具類）

| Skill | Group | 建議 cross-ref hub |
|-------|-------|---------------------|
| `cpp-change` | 程式開發 | scientific-rigor §6 消融原則 |
| `inject-hypothesis` | 研究專用 | scientific-rigor §7.1 Pre-registration |
| `observation-analysis` | 研究專用 | scientific-rigor §5 + auc-confound-guard |
| `results-analysis` | 研究專用 | scientific-rigor §3 Effect Size |
| `review-evidence` | 研究專用 | scientific-rigor §2 Evidence Tier + §8.4 |

## Recommendations

1. **即刻修 (P0)**: 5 YAML invalid → 30 min 機械式 description quote
2. **本週修 (P1)**: 5 個非工具類 D3 cross-ref → 30 min 加 backlink
3. **2 週修 (P2)**: 13 個非工具類 D2 診斷段補完 → 每個 15 min × 13 = 3.5hr
4. **可選 (P3)**: 12 個工具類 D2/D3 → 工具類本就低 cross-ref 需求，pending 評估

## 後續審查階段

- P2 Hooks audit
- P3 Agents audit
- P4 Prompt 庫 / templates audit
- P5 Hooks × Skills × Settings 互動 audit
- P6 Documentation drift audit