# State Directory — CLAUDE.md

> Resilient Waterfall harness 的執行時 state 與 audit 紀錄根目錄。
> Per-folder rules（v1.7 §11.3 batch C 採納）— 此檔僅放 state-specific 約束，root CLAUDE.md 仍是全域權威。

## 結構與職責

```
state/
├── active.json              ← 索引：當前 active cycles（≤5）
├── cycles/                  ← 真實 active cycle 工作目錄
│   └── {cycle_id}/
│       ├── state.json       ← cycle 主狀態（schema 驗證）
│       ├── plan.json
│       ├── precheck.json
│       ├── pilot.json
│       ├── generalize.json
│       ├── evaluation.json
│       └── reflection.log
├── retro_cycles/            ← Drill 1 retro fixture（不在 active.json）
├── cycles_archived/         ← P6 完成後搬移（gitignored 由 cycle）
├── invalidation/            ← stale_marks.jsonl + binary_versions.jsonl
└── schemas/                 ← JSON schema 真理；不可隨意改
```

## 不可手改清單

- **`state.json` / `plan.json` 等 cycle artifact**：必須由對應 skill 寫入（cycle-init / research-loop / check-staleness / run-evaluator / conclude-research）。手改會破壞 schema validation 與 history audit 鏈。
- **`active.json`**：只能由 cycle-init 或 cycle 終結時自動更新；schema 限 maxItems=5。
- **`schemas/*.schema.json`**：改 schema = 改現有 cycle 行為。需 plan v1.x 明確 Decision Log + 不破壞既有 cycle 解析。

## 可手動操作

- **新建 retro_cycles fixture**：手寫 5 件 artifacts 進 retro_cycles/{event_id}/ 用於 regression test（如 Drill 1 retro 模式）。
- **`invalidation/stale_marks.jsonl`** append-only：可手動加 stale 紀錄（hook 也會自動 append）。

## 何時碰 / 不該碰

| 情境 | 應做 | 不該 |
|---|---|---|
| 啟動新 cycle | invoke `/cycle-init` | 直接 mkdir cycles/{id}/ |
| 手動推進 phase | invoke 對應 skill | 直接編輯 state.json `phase` |
| cycle 撤回或 retract | invoke `/conclude-research` 寫 NEGATIVE | 直接刪 cycle 目錄 |
| 想看現狀 | invoke `/cycle-state` | 自己 cat 多個 state.json |

## 與 evidence_ledger 同步規則

`research/autoresearch/evidence_ledger.jsonl` 是**歷史 SoT**，state.json 是**現況快照**。崩潰恢復時可由 ledger 重建 state。**先寫 ledger 再更新 state**（heads-up #4 in plan §10）。

## E1 paths 應用（v1.7）

未來 skill 寫入 state/ 時應在 frontmatter 加：
```yaml
paths: ["state/cycles/**/*.json", "state/active.json"]
```
避免被誤觸發在 docs/ 或 src/ 編輯時。
