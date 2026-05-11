<!--
建立時間: 2026-05-09 03:09
目標: 提供 Codex 在 docs/ 目錄內撰寫、搬移與索引研究文件的局部規則
處理範圍: docs/CURRENT_FOCUS.md、docs/experiments、docs/reports、docs/plans、docs/presentations
關聯檔案:
  - ../AGENTS.md
  - ../docs/CLAUDE.md
  - ../.agents/skills/README.md
-->

# Docs Directory - AGENTS.md

> InterSubMod 研究文檔分流 + INDEX SoT 規範。
> 本檔是 `docs/CLAUDE.md` 的 Codex 版本；`CLAUDE.md` 保留給 Claude Code legacy，Codex 以本檔與 repo root `AGENTS.md` 為準。

## 結構與職責

```
docs/
├── CURRENT_FOCUS.md             <- 當前主軸（研究任務必讀）
├── README.md                    <- 文件導航
├── experiments/
│   ├── INDEX.md                 <- 實驗 SoT；所有 in_progress + validated 都需登錄
│   ├── in_progress/{YYYY}/{MM}/ <- 進行中報告（可改）
│   └── validated/{YYYY}/{MM}/   <- 已 P6 收尾報告（不再改，加 banner 才能改）
├── reports/
│   ├── research_landscape/      <- 9 篇主軸大圖景（半穩定）
│   ├── validated/               <- 已驗證跨 cycle 整合報告
│   └── pi_reports/              <- 給教授的整理
├── plans/                       <- 實驗計劃草稿
├── concepts/                    <- 研究概念樹
├── data_specs/                  <- 資料 schema / format 規格
├── presentations/               <- PPTX 與週報草稿
└── references/                  <- 啟動上下文、手冊
```

## 強制規則

### 命名格式

一般 Markdown 採 `{YYYYMMDD}_{描述}_{NN}.md`，例：`20260507_COLO829_KDE_revisit_01.md`。

- 中英混排不在檔名，用底線分隔。
- 分析報告可保留大寫類別 prefix（如 `O11_`、`S5_`）。
- 路徑含中文字串時，shell 或 git 命令需注意引號保護。

### in_progress vs validated 流程

| 階段 | 路徑 | 何時搬 | 限制 |
|---|---|---|---|
| 草稿 / 探索 | `in_progress/` | cycle 進行中 | 自由改 |
| 已 P6 收尾 | `validated/` | conclude-research 流程搬移 | 加 retraction banner 才能再改 |

手動搬移會觸發 `INDEX.md` 同步更新義務；漏更 INDEX 會讓後續 audit 產生 orphan。

### INDEX.md 是 SoT

`docs/experiments/INDEX.md` 列所有實驗、結論與 tier。任何 `.md` 加進 `docs/experiments/` 必須同步加 INDEX entry。`conclude-research`、週報與 provenance audit 都依賴此 SoT。

## 可手動操作 / 不該碰

| 情境 | 應做 | 不該 |
|---|---|---|
| 寫新實驗報告 | 加在 `in_progress/{YYYY}/{MM}/` 並更新 INDEX | 直接寫 `validated/` |
| 修 validated 報告 | 加 retraction banner + commit message 標 RETRACTED | 直接 edit |
| 看主軸現況 | 讀 `docs/CURRENT_FOCUS.md` | 自己憑記憶整合 INDEX |
| 給教授 review | 寫到 `pi_reports/` 並引用 validated 報告 | 直接傳 in_progress 草稿 |

## 路徑前綴規則

列 `.md` 路徑給使用者時必須以 `InterSubMod/` 前綴，例：

```text
GOOD: InterSubMod/docs/experiments/in_progress/2026/05/20260507_xxx.md
BAD:  docs/experiments/in_progress/2026/05/20260507_xxx.md
```

## Codex 套用方式

- Claude Code slash commands 不直接搬移；舊文件中的 `/doc-standards`、`/weekly-report` 在 Codex 端應改為 `$doc-standards`、`$weekly-report` 或自然語言觸發同名 skill。尚未遷移的 skill 不應假裝可用。
- Claude hooks 不直接搬移；文件 gate 應改成 skill checklist、`INDEX.md` 驗證腳本，或在最終回報列出未驗證風險。
- Claude agents 不直接搬移；文件整合可以使用 Codex subagent，但只有在使用者明確要求平行化或任務已授權分工時使用，且各 subtask 的輸出檔案要互不重疊。

## 未來 paths 提示

若 Codex 後續支援 repo-local skill path filter，docs 相關 skill 可使用：

```yaml
paths: ["docs/**/*.md", "**/*.md"]
```

週報或簡報相關 skill 可使用：

```yaml
paths: ["docs/reports/**/*.md", "docs/presentations/**/*.md"]
```
