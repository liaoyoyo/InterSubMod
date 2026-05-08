# Docs Directory — CLAUDE.md

> InterSubMod 研究文檔分流 + INDEX SoT 規範。
> Per-folder rules（v1.7 §11.3 batch C 採納）。Root CLAUDE.md 為全域權威；本檔只補 docs-specific 約束。

## 結構與職責

```
docs/
├── CURRENT_FOCUS.md             ← 當前主軸（必讀；每對話開頭）
├── README.md                     ← 文件導航
├── experiments/
│   ├── INDEX.md                  ← 實驗 SoT；所有 in_progress + validated 都需登錄
│   ├── in_progress/{YYYY}/{MM}/  ← 進行中報告（可改）
│   └── validated/{YYYY}/{MM}/    ← 已 P6 收尾報告（不再改，加 banner 才能改）
├── reports/
│   ├── research_landscape/       ← 9 篇主軸大圖景（半穩定）
│   ├── validated/                ← 已驗證跨 cycle 整合報告
│   └── pi_reports/               ← 給教授的整理
├── plans/                        ← 實驗計劃草稿
├── concepts/                     ← 研究概念樹
├── data_specs/                   ← 資料 schema / format 規格
├── presentations/                ← PPTX 與週報草稿
└── references/                   ← 啟動上下文、手冊
```

## 強制規則

### 命名格式（doc-standards skill 自動套用）

`{YYYYMMDD}_{描述}_{NN}.md` 例：`20260507_COLO829_KDE_revisit_01.md`

中英混排不在檔名（用底線分隔）；分析報告大寫類別（如 `O11_`, `S5_`）；註：路徑中文字串需 git 引號保護。

### in_progress vs validated 流程

| 階段 | 路徑 | 何時搬 | 限制 |
|---|---|---|---|
| 草稿 / 探索 | `in_progress/` | cycle 進行中 | 自由改 |
| 已 P6 收尾 | `validated/` | conclude-research 自動搬 | 加 retraction banner 才能再改 |

**手動搬移觸發 INDEX.md 同步更新義務** — 漏更 INDEX → 後續 audit 跳出 orphan。

### INDEX.md 是 SoT

`experiments/INDEX.md` 列所有實驗 + 結論 + tier。任何 .md 加進 docs/experiments/ 必須**同步**加 INDEX entry。conclude-research / weekly-report 都依賴此 SoT。

## 可手動操作 / 不該碰

| 情境 | 應做 | 不該 |
|---|---|---|
| 寫新實驗報告 | 加在 `in_progress/{YYYY}/{MM}/` 並更 INDEX | 直接寫 `validated/` |
| 修 validated 報告 | 加 retraction banner + commit message 標 RETRACTED | 直接 edit |
| 看主軸現況 | `Read docs/CURRENT_FOCUS.md` | 自己整合 INDEX |
| 給教授 review | 寫到 `pi_reports/` 並引用 validated 報告 | 直接傳 in_progress 草稿 |

## 路徑前綴規則（與 root CLAUDE.md 對齊）

列 .md 路徑給用戶時必須以 `InterSubMod/` 前綴，例：
```
GOOD: InterSubMod/docs/experiments/in_progress/2026/05/20260507_xxx.md
BAD:  docs/experiments/in_progress/2026/05/20260507_xxx.md
```

## E1 paths 應用（v1.7）

doc-standards skill frontmatter 應加：
```yaml
paths: ["docs/**/*.md", "**/*.md"]
```
weekly-report skill：
```yaml
paths: ["docs/reports/**/*.md", "docs/presentations/**/*.md"]
```
