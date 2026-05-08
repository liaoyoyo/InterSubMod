# Research Directory — CLAUDE.md

> InterSubMod 多週級研究專案 + AutoResearch queue + evidence_ledger SoT。
> Per-folder rules（v1.7 §11.3 batch C 採納）。Root CLAUDE.md 為全域；本檔只補 research-specific 約束。

## 結構與職責

```
research/
├── autoresearch/
│   ├── research_direction.md         ← AutoResearch 候選 queue（候選不執行）
│   ├── hypothesis_queue.json          ← 假說 SoT（inject-hypothesis 寫）
│   └── evidence_ledger.jsonl          ← 歷史 audit 真理（append-only）
├── {project_name}/                    ← 多週專案 scaffolding (init-research 建)
│   ├── manifest.yaml
│   ├── 00_PLAN.md                     ← 專案總計劃
│   ├── observations/
│   ├── scripts/
│   └── ...
├── data_registry/                     ← 跨 project 共用資料索引
└── external_tools/                    ← 第三方 tool 整合腳本
```

## 強制規則

### evidence_ledger.jsonl 是 SoT，不是 state.json

`evidence_ledger.jsonl` 為**歷史軌跡**（append-only，不可 rewrite）；`state/cycles/{id}/state.json` 為**現況快照**。

- 跨 cycle 結論在 ledger 找
- 單 cycle 當下狀態在 state.json 找
- 兩者衝突時：**ledger 為準**（state.json 可由 ledger 重建）

### research/{project} vs state/cycles 分工

| 用 research/ | 用 state/cycles/ |
|---|---|
| 多週長期專案（≥3 假說） | 單假說短驗證（≤2 day） |
| init-research 建立 | cycle-init 建立 |
| 含 manifest.yaml + 多 step | 含 7-phase state machine |
| 範例：v5_provenance_followup | 範例：D2A-COLO829-KDE-rerun |

### research_direction.md 是 queue，不是執行觸發

只放候選方向 + 優先序；**用戶決定哪個動再 inject-hypothesis + cycle-init**。AI 不可自動從 queue 拉項目跑。

### autoresearch 修改規則

- `hypothesis_queue.json`：只能由 inject-hypothesis / pivot-direction skill 寫；手改要 commit message 標 `[manual-queue-edit]`
- `evidence_ledger.jsonl`：**永遠 append-only**；發現舊 entry 錯誤 → 加新 entry 標 `corrects: <old_id>`，不可改舊
- `research_direction.md`：手寫 OK，但加優先序變動需週報層級理由

## 可手動操作 / 不該碰

| 情境 | 應做 | 不該 |
|---|---|---|
| 提新假說 | invoke `/inject-hypothesis` | 直接 edit hypothesis_queue.json |
| 新長期專案 | invoke `/init-research` 建 manifest | 手 mkdir research/foo/ |
| 結論失敗 | invoke `/conclude-research` 寫 NEGATIVE | 直接刪 research/foo/ 或改 ledger entry |
| 跨 cycle 整合 | invoke `/provenance-tier-audit` | 自己寫 query script |

## E1 paths 應用（v1.7）

未來 init-research / inject-hypothesis 等 skill frontmatter 應加：
```yaml
paths: ["research/**/*.md", "research/autoresearch/*.json", "research/autoresearch/*.jsonl"]
```

## 不可 git rm 清單

- `evidence_ledger.jsonl` — append-only 歷史
- `hypothesis_queue.json` — 假說 SoT
- 已驗證 project（如 v5_provenance_followup）的 results

刪除這些屬 Hard Gate（plan §3 + confirmation-protocol）。
