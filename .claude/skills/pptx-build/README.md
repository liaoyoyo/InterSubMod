# pptx-build skill

> InterSubMod PPT 製作 sub-skill（從原 myPPT P1-P5 拆出）。

## 解決什麼

從 weekly-report master_draft.md 或從零開始，產 PPTX + speaker script + audit。
3 層確認 outline → section → slide + 6 模板識別 + §20 主軸聚焦過濾 + Tier 三層分流 + Claude Vision 10-check。

## 何時觸發

| 觸發 | 範例 |
|------|------|
| 簡報 / PPT / PPTX | 「我要做簡報」「教授報告」「投影片」 |
| 從母稿銜接 | weekly-report C4 後 handoff A → `/pptx-build --from-draft <path>` |

## 5 階段流程

```
P1 Audience & Goal       → C1 main thesis ≤ 30 字 ★ fast-track 必停
P2 Outline               → C2 5-7 段
P3 Section batch (×N)    → C3 5-7 slide + focal point ≤ 20 字
P4 Slide build (×24)     → C4 三層分流 + Vision 10-check
P5 Speaker script        → C5 字數 ↔ 時長 ★ fast-track 必停
```

## 與全域 docs 連結

| 主題 | 文件 |
|------|------|
| 5 階段詳細方法論 | `playbook.md` §1-§24 |
| 6 模板 narrative skeleton | `templates/{improvement,comparison,executive,data,concept,academic}_*.md` |
| 5 個 checkpoint prompt | `prompts/*.md` |
| Style library 14 物件 + 12 版型 | `style_library/{objects,layouts,colors,examples}/` |
| v3 case audit 示範 | `audits/v3_self_phasing_24slide_audit.md` |
| Reusable Python tools | `InterSubMod/tools/ppt_toolkit/` |
| Vision 10-check 模板 | `prompts/visual_review.md` |
| --from-draft 介面 | `prompts/from_draft_loader.md` |

## 與其他 skill 關聯

- **上游觸發** ← `weekly-report` C4 handoff A 自動觸發 `/pptx-build --from-draft`
- **平行不重疊** ← `myPPT`（總入口場景識別 → 委派）
- **規範來源** ← `confirmation-protocol` / `doc-standards`
