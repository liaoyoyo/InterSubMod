<!--
建立時間: 2026-05-08 19:35
目標: 記錄 InterSubMod repo-local Codex skills 的來源、範圍與維護規則
處理範圍: .agents/skills 第一批 Claude skills 遷移與後續 Wave 規劃
關聯檔案:
  - ../../AGENTS.md
  - ../../.claude/skills/
-->

# InterSubMod Codex Skills

本目錄是 InterSubMod 的 repo-local Codex skill 入口。Codex 會從 repo root 的 `.agents/skills/` 掃描 `SKILL.md`；`.claude/skills/` 保留為 Claude Code legacy source，不在此遷移中搬移或刪除。

## 第一批遷移

第一批只涵蓋研究循環核心，避免一次載入全部 Claude skills 造成 trigger 擁擠。

| Codex skill | Legacy source | 用途 |
|---|---|---|
| `confirmation-protocol` | `.claude/skills/confirmation-protocol` | 確認時機、Hard Gate、FYI/Review 分級 |
| `known-pitfalls` | `.claude/skills/known-pitfalls` | InterSubMod 已知 AI 陷阱 |
| `auc-confound-guard` | `.claude/skills/auc-confound-guard` | AUC confound 三關驗證 |
| `cycle-state` | `.claude/skills/cycle-state` | active cycles 唯讀 dashboard |
| `cycle-init` | `.claude/skills/cycle-init` | P0 REGISTER cycle 初始化 |
| `research-loop` | `.claude/skills/research-loop` | P1 PLAN 研究迴圈設計 |
| `validation-protocol` | `.claude/skills/validation-protocol` | L1-L4 假說驗證層級；因 `research-loop` 相依而補入 |
| `check-staleness` | `.claude/skills/check-staleness` | P2 PRECHECK freshness gate |
| `feature-layered-observation` | `.claude/skills/feature-layered-observation` | P3 PILOT feature observation |
| `multi-sample-consistency` | `.claude/skills/multi-sample-consistency` | P4 cross-sample consistency |
| `run-evaluator` | `.claude/skills/run-evaluator` | P5 EVALUATE retraction risk |
| `results-report` | `.claude/skills/results-report` | 實驗結果報告 |
| `conclude-research` | `.claude/skills/conclude-research` | P6 COMMIT / 研究收尾 |

## Codex 維護規則

- `SKILL.md` frontmatter 只保留 `name` 與 `description`，讓 Codex implicit trigger 能穩定匹配。
- 若內容提到 `/skill-name`，在 Codex 中等同 `$skill-name` 或同名 skill 明確觸發；保留 slash 寫法只為相容既有研究文件。
- 支援檔案只搬第一批實際引用的 `scripts/`、`references/`、`examples/` 與資料表；大型 PPT、weekly report、citation resources 不在第一批。
- 不直接刪除或覆寫 `.claude/skills`。若 legacy source 後續更新，要人工比較後再同步到 `.agents/skills`。
- 網路政策：本地 repo docs、Knowledge Base、MCP 優先；web search 只作 fallback；live network 需要限 domain 並標註來源。

## Wave 2 候選

下一批可評估：

- `weekly-report`
- `pptx-build`
- `myPPT`
- `citation-verification`
- `problem-framing-ideation`
- `structured-tech-report`

Wave 2 需要額外處理互動式確認、PPT 資源、citation web verification、以及可能的 plugin 化。
