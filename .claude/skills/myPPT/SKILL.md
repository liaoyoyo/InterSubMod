---
name: myPPT
description: InterSubMod 整體 PPT 工作流總入口（從資料搜集問答到產出可顯示 PPT）。輕量場景識別 → 委派 sub-skill：週報情境委派 weekly-report (W1-W7 母稿)，已有母稿/從零開始委派 pptx-build (P1-P5 PPTX 製作)。本 skill 不直接做母稿/PPTX，只做場景判斷與路徑分流。觸發：「我要做簡報」「教授報告」「週報」「PPT」「向教授彙報」「PPTX」「報告」「整理本週」「投影片」
allowed-tools: Read, AskUserQuestion
user-invocable: true
---

# myPPT — 整體 PPT 工作流總入口

從用戶意圖識別場景，委派對應 sub-skill：

```
[U: 我要做簡報/週報/教授報告]
        ↓
   myPPT 場景識別
        ↓
   AskUserQuestion 4 options
        ↓
   ┌────┬────┬────┬────┐
   ↓    ↓    ↓    ↓
  週報   有母稿  從零  混合（先週報後 PPT）
   ↓    ↓    ↓    ↓
  /weekly-report  /pptx-build  /pptx-build  weekly-report → handoff A → pptx-build
                  --from-draft
```

## 場景識別 4 options（核心邏輯）

```yaml
question: "您要做什麼類型的簡報/報告？"
header: "場景識別"
options:
  - label: "週報整理（給教授/PI 看的進度）"
    description: "→ 委派 weekly-report skill 走 W1-W7 母稿。完成後可選 handoff A 接 pptx-build 產 PPTX"
  - label: "已有母稿（master_draft.md）→ 直接產 PPT"
    description: "→ 委派 pptx-build --from-draft <path>，自動讀 frontmatter 跳過 P1 部分項目"
  - label: "從零開始做簡報（非週報情境）"
    description: "→ 委派 pptx-build 走完整 P1-P5（教授級報告 / lab meeting / 學術 defense）"
  - label: "完整 pipeline（先週報母稿，再產 PPT）"
    description: "→ 先 weekly-report 走 W1-W7，C4 後自動 handoff A 觸發 pptx-build"
```

## 委派邏輯

| 場景 | 委派 | 觸發指令 |
|------|------|---------|
| 週報整理 | weekly-report | （直接觸發 weekly-report skill）|
| 已有母稿 | pptx-build | `/pptx-build --from-draft <path>` |
| 從零做簡報 | pptx-build | （直接觸發 pptx-build skill）|
| 完整 pipeline | weekly-report → pptx-build | weekly-report C4 handoff A 自動接 pptx-build |

## 與 sub-skill 關聯

- **weekly-report** — W1-W7 母稿生成（4 主線類型 / 4 層分類 [F]/[O]/[I]/[U] / 教授追問預測 / Layer 0-4 + 17 段 + §0 Highlights）
- **pptx-build** — P1-P5 PPTX 製作（6 模板識別 / §20 主軸聚焦 6 階段過濾 / Tier 三層分流 / Vision 10-check / style_library + ppt_toolkit）

## 不直接做

- ❌ raw data 收集（→ weekly-report W1）
- ❌ 母稿撰寫（→ weekly-report W7）
- ❌ outline / section / slide build（→ pptx-build P2-P4）
- ❌ Vision 10-check（→ pptx-build P4）

## 詳見

- 整體 pipeline 設計 → `playbook.md`
- workflow 圖 → `workflow_diagrams.md`
- 三 skill 角色分工 → `README.md`
- 母稿 schema → `InterSubMod/.claude/skills/weekly-report/references/HANDOFF_TO_PPTX_BUILD.md`

## 注意事項

1. 不要直接在 myPPT skill 內做實作 — 必須委派
2. 觸發詞與 weekly-report / pptx-build 部分重疊，由 myPPT 優先觸發場景識別
3. 場景識別 AskUserQuestion 為**必停**（不論 fast-track）
