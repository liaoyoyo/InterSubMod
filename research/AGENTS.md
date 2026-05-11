<!--
建立時間: 2026-05-09 03:09
目標: 提供 Codex 在 research/ 目錄內操作多週研究專案與 AutoResearch ledger 的局部規則
處理範圍: research/autoresearch、research/{project_name}、research/data_registry、research/external_tools
關聯檔案:
  - ../AGENTS.md
  - ../research/CLAUDE.md
  - ../.agents/skills/README.md
-->

# Research Directory - AGENTS.md

> InterSubMod 多週級研究專案 + AutoResearch queue + evidence_ledger SoT。
> 本檔是 `research/CLAUDE.md` 的 Codex 版本；`CLAUDE.md` 保留給 Claude Code legacy，Codex 以本檔與 repo root `AGENTS.md` 為準。

## 結構與職責

```
research/
├── autoresearch/
│   ├── research_direction.md         <- AutoResearch 候選 queue（候選不執行）
│   ├── hypothesis_queue.json         <- 假說 SoT（inject-hypothesis 寫）
│   └── evidence_ledger.jsonl         <- 歷史 audit 真理（append-only）
├── {project_name}/                   <- 多週專案 scaffolding（init-research 建）
│   ├── manifest.yaml
│   ├── 00_PLAN.md                    <- 專案總計劃
│   ├── observations/
│   ├── scripts/
│   └── ...
├── data_registry/                    <- 跨 project 共用資料索引
└── external_tools/                   <- 第三方 tool 整合腳本
```

## 強制規則

### evidence_ledger.jsonl 是 SoT，不是 state.json

`evidence_ledger.jsonl` 是歷史軌跡，append-only，不可 rewrite；`state/cycles/{id}/state.json` 是現況快照。

- 跨 cycle 結論在 ledger 找。
- 單 cycle 當下狀態在 state.json 找。
- 兩者衝突時：ledger 為準，state.json 可由 ledger 重建。

### research/{project} vs state/cycles 分工

| 用 research/ | 用 state/cycles/ |
|---|---|
| 多週長期專案（>=3 假說） | 單假說短驗證（<=2 day） |
| `init-research` 建立 | `cycle-init` 建立 |
| 含 `manifest.yaml` + 多 step | 含 7-phase state machine |
| 範例：`v5_provenance_followup` | 範例：`D2A-COLO829-KDE-rerun` |

### research_direction.md 是 queue，不是執行觸發

只放候選方向與優先序；使用者決定哪個方向要啟動後，才使用 `inject-hypothesis` 與 `cycle-init`。AI 不可自動從 queue 拉項目執行。

### autoresearch 修改規則

- `hypothesis_queue.json`：只能由 `inject-hypothesis` / `pivot-direction` skill 寫；手改要 commit message 標 `[manual-queue-edit]`。
- `evidence_ledger.jsonl`：永遠 append-only；發現舊 entry 錯誤時，加新 entry 標 `corrects: <old_id>`，不可改舊。
- `research_direction.md`：手寫 OK，但優先序變動需有週報層級理由。

## 可手動操作 / 不該碰

| 情境 | 應做 | 不該 |
|---|---|---|
| 提新假說 | 使用 `$inject-hypothesis` 或同名 Codex skill | 直接 edit `hypothesis_queue.json` |
| 新長期專案 | 使用 `$init-research` 建 manifest | 手 mkdir `research/foo/` |
| 結論失敗 | 使用 `$conclude-research` 寫 NEGATIVE | 直接刪 `research/foo/` 或改 ledger entry |
| 跨 cycle 整合 | 使用 `$provenance-tier-audit` | 自己寫一次性 query script |

## 不可 git rm 清單

- `evidence_ledger.jsonl`：append-only 歷史。
- `hypothesis_queue.json`：假說 SoT。
- 已驗證 project（如 `v5_provenance_followup`）的 results。

刪除這些屬 Hard Gate；若需要移除，依 repo `AGENTS.md` 的 archive 規則，先確認再搬移。

## Codex 套用方式

- Claude Code slash commands 不直接搬移；舊文件中的 `/inject-hypothesis`、`/init-research`、`/provenance-tier-audit` 在 Codex 端應改為 `$inject-hypothesis`、`$init-research`、`$provenance-tier-audit` 或自然語言觸發同名 skill。
- Claude hooks 不直接搬移；ledger、queue、state 的 gate 應改成 append-only 規則、schema 驗證腳本與 skill 內 checklist。
- Claude agents 不直接搬移；長期研究拆工時只在使用者授權下使用 Codex subagent，且每個 subtask 要有清楚 ownership 與互不重疊的寫入範圍。

## 未來 paths 提示

若 Codex 後續支援 repo-local skill path filter，research 相關 skill 可使用：

```yaml
paths: ["research/**/*.md", "research/autoresearch/*.json", "research/autoresearch/*.jsonl"]
```
