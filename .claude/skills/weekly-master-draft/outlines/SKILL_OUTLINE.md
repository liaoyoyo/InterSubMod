---
title: weekly-master-draft SKILL.md 大綱（待用戶審查）
date: 2026-05-01
status: outline-for-review
target_lines: 150-180
---

# weekly-master-draft SKILL.md 大綱

## 核心定位

**「不是做 PPT。是做 PPT 之前的『資料佐證母稿』。」**

母稿 = 經過 AI 整理 + 用戶 4 層分類確認 + 邏輯檢查 + 教授視角預演的 17 段結構化 .md 文件。
完成後 → 自動建議銜接 myPPT 把母稿做成 deck（用戶可選擇繼續或停在母稿）。

## 9 個區塊規劃

### 區塊 1：YAML frontmatter（≈ 8 行）

```yaml
---
name: weekly-master-draft
description: Use whenever the user wants to organize weekly research progress for advisor reporting (週報/週進度/向教授報告/PI 週彙報/weekly progress/research weekly). Builds a 17-section structured master draft (.md) BEFORE any PPT creation, by separating fact/observation/inference/unconfirmed, identifying the report's main thread (progress/problem/consult/exploration), running 5-7 professor-question prediction, and checking overclaim/laundry-list red flags. Output is a data-backed master draft .md that can be handed off to myPPT skill for slide generation.
allowed-tools: Read, Write, Edit, AskUserQuestion, Glob, Grep
---
```

**觸發 keyword（強度 ≥ 9/10）**：週報、weekly report、週進度、向教授報告、PI 週彙報、研究進度、本週進展、weekly meeting

### 區塊 2：MUST-INVOKE-FIRST 警告（≈ 12 行）

> 本 skill 的核心原則：**先有母稿，才有簡報。**
>
> 強制執行：
> 1. 不可跳過內容 4 層分類（事實/觀察/推論/待確認）
> 2. 不可一次問用戶 >5 個問題（每輪 3-5 題）
> 3. 不可把推論寫成事實
> 4. 不可省略教授問答預測
> 5. 不可直接跳到產 PPT（先完成母稿，再用 AskUserQuestion 詢問是否銜接 myPPT）

### 區塊 3：7 階段流程速覽（≈ 25 行）

```
[U: 我要做週報]
    ↓
W1. Raw Data 收集            → C0 確認原始進度資料來源
    ↓
W2. 主線類型識別             → C1 進展型/問題型/求協助型/探索型
    ↓
W3. 內容 4 層分類            → C2 事實/觀察/推論/待確認 標籤
    ↓
W4. 重點排序 + 分流          → C3 PPT/講稿/備註/暫存 4 桶
    ↓
W5. 邏輯檢查（紅旗）          → C4 流水帳/過度宣稱/教授視角缺
    ↓
W6. 教授問答預測             → C5 5-7 個追問 + 預備回答
    ↓
W7. 17 段母稿產出            → C6 母稿最終確認
    ↓
[Output: weekly_master_draft_YYYYMMDD.md]
    ↓
[Q via AskUserQuestion: 要繼續銜接 myPPT 產 PPT 嗎？]
    ├─ 是 → 觸發 myPPT skill，自動讀取母稿
    └─ 否 → 母稿留檔，週報結束
```

### 區塊 4：4 種主線類型識別（≈ 18 行）

| 觸發 keyword/語境 | 主線類型 | 推薦敘事弧 |
|----------------|---------|---------|
| 「進展、突破、完成、達成」| **進展型** | 背景 → 本週處理 → 結果 → 初步分析 → 下週 |
| 「問題、卡住、blocker、bug、anomaly」| **問題型** | 方法 → 問題發現 → 目前判斷 → 求建議 |
| 「不確定、需要 advisor、方向選擇」| **求協助型** | 情境 → 多選項 → 各選項利弊 → 待決策點 |
| 「新方向、pilot、初步觀察、探索」| **探索型** | 動機 → 假設 → pilot 結果 → 是否值得投入 |

混合用例規則：以「教授最關心的點」為主軸，其他敘事弧降為 sub-thread。

### 區塊 5：強制 checkpoint（C0-C6）清單（≈ 30 行）

7 個必停點：
- **C0 Raw Data 確認**：本週實際做了什麼（時間軸 + 產出列表）
- **C1 主線類型**：4 選 1 + 一句話 main statement（≤ 30 字）
- **C2 內容 4 層分類**：每筆素材標 [F]/[O]/[I]/[U]
- **C3 重點排序**：對每筆素材決定 PPT / 講稿 / 備註 / 暫存
- **C4 邏輯紅旗檢查**：通過/失敗清單，失敗項回到 C2 重分類
- **C5 教授問答預演**：5-7 個追問 + 預備回答
- **C6 17 段母稿確認**：最終 .md 用戶批准

每個 checkpoint 含：暫停理由 + 預期輸出 + AskUserQuestion 模板（指向 prompts/）

### 區塊 6：紅旗清單（≈ 25 行）

**過度宣稱紅旗**（C4 觸發）：
- 「證實」「確認」「解決」用於初步觀察 → 改「初步觀察、需 N 樣本驗證」
- 「全部」「完全」用於部分樣本結果 → 改具體 sample 範圍
- 「顯著」未含 p-value 或 CI → 補統計或改「具方向性」
- 主動語態斷言而無證據佐證 → 補來源或降為推論

**流水帳紅旗**（C4 觸發）：
- 「本週做了 ABCDEFG」（>5 件平列）→ 重排優先序，>3 個降到備註
- 無因果連接詞（因為/所以/由於）→ 改寫成因果鏈
- 每件事獨立段落無串接 → 至少 2 件合併為一個 narrative

**教授視角缺紅旗**（C5 觸發）：
- 母稿無「教授可能問」段 → 強制補
- 母稿無「下週計畫銜接本週發現」 → 強制補
- 母稿無「需要教授判斷的點」（求協助型必要）→ 強制補

### 區塊 7：與 myPPT 的銜接（≈ 15 行）

**輸出協議（contract）**：
- 母稿輸出路徑：`InterSubMod/docs/weekly_reports/YYYYMMDD/master_draft_v{N}.md`
- 母稿 frontmatter 含：`report_type`、`main_thread`、`audience`、`target_duration_min`、`source_artifacts: [...]`
- myPPT 讀取母稿後，跳過 §1 Audience & Goal 的「main thesis 鎖定」（已由 C1 完成）

**handoff prompt（C6 後）**：
> 母稿已完成（17 段，N 字）。是否要繼續銜接 myPPT 產 PPT？
> - 是 → 觸發 myPPT，自動帶入母稿路徑
> - 否 → 母稿留檔，週報任務結束
> - 是但要先休息 → 母稿留檔，下次再用 `/myPPT --from-draft <path>` 啟動

### 區塊 8：fast-track 例外（≈ 12 行）

用戶明示「全自動」/「auto」時：
- 跳過 C2 / C3 逐項確認（保留 C0 / C1 / C5 / C6）
- 4 層分類用 AI 預設標籤（用戶最終 review）
- 教授問答預測仍強制（不可省）

但 Hard Gate 仍維持：
- 主線類型必選 4 之一（不可空白）
- 「事實」標籤需有具體來源檔案路徑（無 source = 強制降「觀察」）

### 區塊 9：與 myPPT / 其他 skill 的關聯（≈ 12 行）

| Skill | 觸發時機 | 關係 |
|-------|---------|------|
| myPPT | C6 後用戶選擇繼續 | **下游接棒** |
| weekly-report（既有） | 用戶要求「週報但不報告」 | 重疊：weekly-report 偏向 evidence ledger 整理；weekly-master-draft 偏向「向教授彙報的論述」。可考慮 weekly-report 升級為 weekly-master-draft 的 W1 raw data 收集器 |
| structured-tech-report | 單一工程改動 deep dive | **平行**（不同 cadence） |
| confirmation-protocol | C0-C6 對應 Gate 級別 | **規範來源** |

## 大綱審查重點

1. 9 區塊組合是否合理？
2. 觸發 keyword 是否完整？是否需要加「lab meeting」「dry run」「review meeting」？
3. 7 階段（W1-W7）是否清楚？要 ASCII 還是 mermaid？
4. C0-C6 七個 checkpoint 是否過多？（vs myPPT 的 C1-C5 五個）能否合併？
5. 「過度宣稱」紅旗是否完整？是否需要再列 4-5 條？
6. handoff 路徑 `InterSubMod/docs/weekly_reports/YYYYMMDD/` 是否合適？還是用 `docs/presentations/weekly/`？
7. fast-track 條件是否符合您的「全自動」習慣？
8. 與既有 weekly-report skill 是否衝突？要合併還是並行？
