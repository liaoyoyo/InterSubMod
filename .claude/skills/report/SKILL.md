---
name: report
description: AI 對話執行報告撰寫。記錄重要對話的決策、程式碼修改、後續行動。USE WHEN：會話結束前撰寫執行報告、「寫報告」「session report」「conversation log」、Stop hook 提醒時。輸出到 docs/provenance/ai_sessions/*.md。SKIP WHEN 寫週報（用 weekly-report）、寫實驗結果報告（用 results-report）、寫 13 段技術報告（用 structured-tech-report）、寫研究收尾結論（用 conclude-research）、純 build / commit / docs 寫作。
allowed-tools: Read, Write, Glob, Grep
user-invocable: true
---

> **⚠ 2026-05-22 thin wrapper migration**: 本 skill 為 `/narrative-frame` skill 的 **thin wrapper**。
> 預設 framework = **AI-Session-Companion**（Timeline + Key decisions + Provenance）。
> 等同 `/narrative-frame apply AI-Session-Companion`。
> 用戶可隨時換 framework：對話中説「用 13 段技術報告風」或直接走 `/narrative-frame N1-N6` 動態挑。
> Catalog: `InterSubMod/.claude/skills/narrative-frame/references/framework_catalog.md` §11 hybrid。

# AI 對話報告撰寫技能

你是一位專業的技術文件撰寫者，負責記錄 AI 對話的重要內容。

## 使用時機

在以下情況使用此技能：
- 完成重要的開發任務後
- 做出關鍵的技術決策後
- 會話即將結束前
- 需要記錄重要的觀念改變

## 報告內容

1. **對話目標**：本次對話要完成的任務
2. **關鍵決策**：過程中做出的重要決定
3. **修改清單**：所有被修改的檔案
4. **觀念更新**：任何概念或想法的改變
5. **後續行動**：建議的下一步

## 輸出位置

`docs/provenance/ai_sessions/{YYYY}/{MM}/`

## 檔案命名

`{YYYYMMDD}_{對話主題}_執行報告_01.md`

## 報告模板

```markdown
# {對話主題} 執行報告

<!--
建立時間: YYYY-MM-DD HH:MM
目標: [對話目標]
處理範圍: [涵蓋的工作範圍]
關聯檔案:
  - [相關檔案路徑 1]
  - [相關檔案路徑 2]
-->

## 對話資訊

| 項目 | 內容 |
|------|------|
| 日期 | YYYY-MM-DD |
| 主要目的 | ... |
| AI 模型 | Claude Opus 4.5 |

## 對話背景

[簡要描述對話的背景和起因]

## 關鍵決策

| 決策 | 原因 | 影響 |
|------|------|------|
| ... | ... | ... |

## 產出成果

### 新增檔案
- `path/to/file` - 說明

### 修改檔案
- `path/to/file` - 說明

### 刪除檔案
- `path/to/file` - 說明

## 觀念更新

- [更新前] → [更新後]

## 後續行動

- [ ] 行動 1
- [ ] 行動 2

## 對話摘要

[用 3-5 句話總結本次對話的主要內容和成果]
```

## 撰寫指南

1. **保持客觀**：如實記錄，不加入主觀評價
2. **重點突出**：聚焦在重要的決策和變更
3. **可追溯**：確保記錄的檔案路徑準確
4. **可執行**：後續行動應該具體可執行
5. **簡潔明瞭**：避免冗長的描述
