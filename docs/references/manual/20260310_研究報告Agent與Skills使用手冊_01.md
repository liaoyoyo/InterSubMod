<!--
建立時間: 2026-03-10 12:55
目標: 說明如何使用 InterSubMod 研究文件 Agent 與三個 Skills，並提供其他人可直接使用的指令模板
處理範圍: 週報、證據鏈、背景整理、偏好修正
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/.claude/agents/intersubmod-weekly-research-agent.md
  - /home/liaoyoyo2001/.codex/skills/intersubmod-context-synthesizer/SKILL.md
  - /home/liaoyoyo2001/.codex/skills/intersubmod-weekly-report-writer/SKILL.md
  - /home/liaoyoyo2001/.codex/skills/intersubmod-report-prompt-refiner/SKILL.md
-->

# 研究報告 Agent 與 Skills 使用手冊

## 1. 快速開始

### 建議使用順序

1. 讀 [CURRENT_FOCUS.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md)
2. 讀 [INDEX.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/INDEX.md)
3. 讀 [README.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/README.md)
4. 交給 [intersubmod-weekly-research-agent.md](/big8_disk/liaoyoyo2001/InterSubMod/.claude/agents/intersubmod-weekly-research-agent.md) 或直接指定對應 skill

## 2. 什麼情況用哪個技能

| 任務 | 建議工具 |
|---|---|
| 補背景、寫導讀、術語整理 | `intersubmod-context-synthesizer` |
| 產出週報或證據鏈主報告 | `intersubmod-weekly-report-writer` |
| 修正 prompt、收斂偏好 | `intersubmod-report-prompt-refiner` |
| 依廖子游固定風格做研究簡報 / PPTX | `liao-research-ppt` + `pptx` |
| 一次要做完整研究文件整理 | `intersubmod-weekly-research-agent` |

## 3. 目前仍需補齊的內容

### 研究內容缺口

1. 是否需要固定的對外摘要版格式。
2. 是否要針對圖表建立更強制的 caption 模板。

### 證據鏈缺口

1. 某些跨樣本結論仍可能缺正式 validated report。
2. 某些 diagnostics 圖片尚未全部配好文字解釋。

### 指令規格缺口

1. 是否要固定輸出篇幅上限。
2. 是否要區分內部深度版與管理層精簡版。
3. 是否要固定簡報輸出到 `docs/presentations/` 並附版本確認檔。

## 3.1 若任務是做 PPTX

先讀：

1. [個人 PPT 設計風格規範](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260311_個人PPT設計風格規範_01.md)
2. [個人 PPT profile 設定檔](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/assets/20260311_liao_research_ppt_profile_01.json)
3. [liao-research-ppt skill](/home/liaoyoyo2001/.codex/skills/liao-research-ppt/SKILL.md)
4. [研究週報 PPTX 客製化設定與製作手冊](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260311_研究週報PPTX客製化設定與製作手冊_01.md)

輸出規則：

1. `.md` 報告留在 `docs/reports/...`
2. `.pptx` 以 `docs/presentations/...` 作為檢視入口
3. 每份 deck 要附版本確認檔

## 4. 其他人如何撰寫指令

### 簡短版

適用：已經很熟悉脈絡，只要快速整理。

```text
請使用 intersubmod-weekly-report-writer，整理 2026-03-05 到 2026-03-10 的研究週報，先講重點結論，再附 benchmark、規則比較與下一步。
```

### 完整版

適用：希望有背景、時間序與主題整合。

```text
請先閱讀 /big8_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md、
/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/INDEX.md、
/big8_disk/liaoyoyo2001/InterSubMod/docs/README.md，
再使用 intersubmod-context-synthesizer 與 intersubmod-weekly-report-writer，
整理指定期間的研究週報。請用台灣繁體中文，開頭先寫重點結論，
所有 benchmark 要列 truth_total、TP、FP、FN、precision、recall、F1，
所有重要檔案請用絕對路徑 Markdown 連結。
```

### 高精度版

適用：需要證據鏈、表格規則、待補齊事項與下次 prompt。

```text
請使用 intersubmod-weekly-research-agent。
範圍：2026-03-05 到 2026-03-10。
主題：TO rescue、methylation-support、caller-first、artifact-veto。
讀者：已理解知識庫內容的研究同仁。
輸出：混合版研究週報。
要求：
1. 開頭先用一段高密度結論破題。
2. 補必要背景知識。
3. 依時間順序整理進展，再做主題式整合分析。
4. 用表格清楚比較 benchmark 與 rescue 規則。
5. 圖片若引用，必須解釋 X/Y 軸意義。
6. 結尾列待補齊事項、下一步與下次可重用 prompt。
```

## 5. 問答式修正範本

若你要持續修正 agent/skill，建議直接回答以下問題：

1. 這次背景知識是太多、太少，還是剛好？
2. 你要更多數據表，還是更多推論文字？
3. 你要的讀者是研究內部版，還是對外摘要版？
4. 哪些表格欄位是一定要固定保留的？
5. 哪些偏好應該升級成固定規則？
6. 這份簡報要偏研究週報模板還是研究驗證模板？
7. 是否需要完整 speaker notes？
8. 是否要產生 deck 版本確認檔與附錄路徑頁？

## 6. Skill / Agent 草案如何繼續演進

下一輪若要更像正式 skill/agent，可優先補：

1. 真實使用案例與對照 prompt
2. 固定模板檔
3. 小型 evaluation prompts
4. prompt 修正 changelog

## 7. 建議下一步

1. 先用這套 agent/skills 寫一份本週週報。
2. 依你的回饋補強固定 prompt。
3. 再決定要不要把 prompt refinement 變成正式 changelog 或 eval 流程。
