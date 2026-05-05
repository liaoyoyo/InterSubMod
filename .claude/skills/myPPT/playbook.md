# myPPT playbook — 整體 pipeline 主方法論

> 本檔說明三 skill 串接邏輯與場景判斷規則。
> 不重複 weekly-report / pptx-build sub-skill 細節（請至各 sub-skill 內 playbook 查閱）。

## §1 整體 pipeline 概念

從用戶意圖到 PPTX 產出的完整路徑：

```
意圖：「我要做簡報/週報/教授報告」
    ↓
階段 0：場景識別（myPPT 本 skill）
    ↓ 委派
階段 W：母稿生成（weekly-report sub-skill，W1-W7）
    ↓ handoff
階段 P：PPTX 製作（pptx-build sub-skill，P1-P5）
    ↓
產出：PPTX + speaker script + audit
```

**核心原則**：先有母稿，再有 PPT。

## §2 場景識別 4 種

### 場景 A：週報整理（最常見）

- 觸發：「整理本週」「週報」「向教授彙報」「PI 月會」
- 委派：weekly-report
- 後續：用戶 C4 後 handoff 4 選決定是否接 pptx-build

### 場景 B：已有母稿，產 PPT

- 觸發：用戶提供 master_draft.md path + 「產 PPT」
- 委派：pptx-build --from-draft <path>
- 跳過 P1 main thesis 鎖定（母稿 frontmatter 已有）
- 跳過 P1 報告類型識別（從 report_type 自動載入）

### 場景 C：從零做簡報（非週報）

- 觸發：「教授級報告」「lab meeting」「academic defense」
- 委派：pptx-build 走完整 P1-P5
- 自動載入 academic_defense / executive_summary 等模板

### 場景 D：完整 pipeline（一條龍）

- 觸發：「我要做週報，最後也要 PPT」
- 委派序列：weekly-report → C4 handoff A → pptx-build --from-draft
- 兩 sub-skill 連續跑，中間用戶可在 C4 暫停

## §3 場景識別決策樹

```
用戶意圖
   ├─ 含「週報」/「整理本週」/「向教授」？
   │   ├─ 是 → 場景 A 或 D
   │   │   └─ 問用戶：「最後要產 PPTX 嗎？」
   │   │       ├─ 是 → 場景 D（完整 pipeline）
   │   │       └─ 否 → 場景 A（只產母稿）
   │   └─ 否 → 進下一條件
   │
   ├─ 提供 master_draft.md path？
   │   ├─ 是 → 場景 B（已有母稿）
   │   └─ 否 → 進下一條件
   │
   ├─ 含「PPT」/「投影片」/「deck」/「教授報告」？
   │   ├─ 是 → 場景 C（從零做簡報）
   │   └─ 否 → AskUserQuestion 4 選
```

## §4 委派指令對應

| 場景 | 委派指令 / 流程 |
|------|---------------|
| A | 直接觸發 weekly-report skill（不需指令）|
| B | `/pptx-build --from-draft <path>` |
| C | 直接觸發 pptx-build skill |
| D | weekly-report → C4 handoff A → 自動 `/pptx-build --from-draft` |

## §5 場景識別 AskUserQuestion 標準模板

詳見 `SKILL.md` 中「場景識別 4 options」段。**必停**（不論 fast-track）。

理由：場景判斷錯誤 = 後續整條 pipeline 走錯，不可省。

## §6 不重複 sub-skill 內容

myPPT 本 skill **只負責場景識別與委派**，不重複：

- W1-W7 細則 → `InterSubMod/.claude/skills/weekly-report/playbook.md`
- W7 母稿格式 + §0 Highlights → `InterSubMod/.claude/skills/weekly-report/references/LAYER_STRUCTURE.md`
- handoff 4 選 → `InterSubMod/.claude/skills/weekly-report/references/HANDOFF_TO_PPTX_BUILD.md`
- P1-P5 細則 → `InterSubMod/.claude/skills/pptx-build/playbook.md`
- §20 主軸聚焦 6 階段 → `InterSubMod/.claude/skills/pptx-build/playbook.md` §20
- 6 報告模板 → `InterSubMod/.claude/skills/pptx-build/templates/`
- style_library + ppt_toolkit → `InterSubMod/.claude/skills/pptx-build/style_library/` + `InterSubMod/tools/ppt_toolkit/`

## §7 與其他 skill 互動

- **structured-tech-report**：平行（單一工程改動 deep dive，非週期）
- **confirmation-protocol**：規範來源（場景識別必停 vs sub-skill checkpoint）
- **doc-standards**：母稿 / PPTX 命名規範

## §8 Failure modes

| 情境 | 處置 |
|------|------|
| 用戶意圖模糊（無明確觸發詞） | AskUserQuestion 4 選讓用戶確認 |
| 用戶選 B 但 master_draft.md 不存在 | 提示「母稿不存在，改場景 C 從零做？」|
| 用戶選 D 但 weekly-report 已 C4 留檔 | 提示「母稿已存於 X，改場景 B 直接接 pptx-build？」|
| 用戶混用觸發詞（如「PPT 週報」）| AskUserQuestion 確認，預設場景 D |

## §9 success criteria

| 指標 | 通過標準 |
|------|---------|
| 場景識別不誤判 | A/B/C/D 場景測試各 1 次，識別正確率 100% |
| 委派完整 | 識別後正確 trigger sub-skill，不直接做實作 |
| 不重複 sub-skill 細節 | 本 skill SKILL.md / playbook 不含 W1-W7 / P1-P5 細則 |
| 觸發優先序 | 「週報」「PPT」混合觸發詞時，myPPT 優先於 sub-skill |
