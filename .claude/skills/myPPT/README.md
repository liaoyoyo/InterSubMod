# myPPT skill (整體 pipeline 總入口)

> InterSubMod PPT 工作流的輕量總入口 — 場景識別 + 委派 sub-skill。

## 三 skill 架構

```
myPPT (總入口，本 skill)
   ├─ 場景識別 → 委派
   │
   ├─ weekly-report (sub-skill: 母稿生成)
   │   W1 raw data → W2 主線 → W3 4 層分類 → W4 排序
   │   → W5 紅旗 → W6 教授追問 → W7 Layer 0-4 + 17 段
   │   → C4 後 handoff 4 選
   │
   └─ pptx-build (sub-skill: PPTX 製作)
       P1 audience → P2 outline → P3 section
       → P4 slide build (Vision 10-check) → P5 speaker script
       → PPTX + speaker + audit
```

## 何時觸發 myPPT

| 觸發語境 | 範例 |
|----------|------|
| 「我要做簡報」/ 不確定怎麼做 | 「我下週要報告」「教授要看簡報」 |
| 「PPT」/ 「投影片」 | 「製作 PPT」「投影片」「deck」 |
| 「週報」/ 「向教授彙報」 | 「整理本週」「PI 月會」 |

myPPT 觸發後必跑場景識別 AskUserQuestion，**不直接做實作**。

## 何時不觸發 myPPT（直接觸發 sub-skill）

- 用戶明確說「跑 weekly-report」/「整理本週進度」 → 直接 weekly-report
- 用戶提供母稿 path + 「產 PPT」 → 直接 pptx-build --from-draft
- 用戶在已有 conversation 中接續 PPT 製作 → 直接 pptx-build

## 與全域 docs 連結

- 場景識別流程 → `SKILL.md`
- 整體 pipeline 設計 → `playbook.md`
- workflow 圖（mermaid）→ `workflow_diagrams.md`
- weekly-report sub-skill → `InterSubMod/.claude/skills/weekly-report/`
- pptx-build sub-skill → `InterSubMod/.claude/skills/pptx-build/`
- 三 skill 架構決策紀錄 → `InterSubMod/.claude/skills/weekly-master-draft/outlines/UPGRADE_PLAN_FOR_WEEKLY_REPORT.md`

## v1 → v2 重大重組（2026-05-05）

舊版 myPPT 是「單一 PPT 製作 skill」，包含 P1-P5 全流程 + style_library。v2 拆分：
- 母稿生成部分 → `weekly-report` sub-skill（升級為 W1-W7）
- PPTX 製作部分 → `pptx-build` sub-skill（保留 P1-P5）
- 場景識別部分 → 本 skill（輕量總入口）

理由：用戶意圖中「先有母稿再有 PPT」是核心邏輯，舊版 myPPT 直接從 outline 開始忽略了「論述母稿」階段。詳見 `InterSubMod/.claude/skills/weekly-master-draft/outlines/UPGRADE_PLAN_FOR_WEEKLY_REPORT.md`。
