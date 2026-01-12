<!--
最後更新: 2026-01-12 00:00
-->

# 當前目標

## 進行中

（目前無進行中的主要目標）

## 已完成（本週）

- 文檔庫架構重整（第二階段）- 2026-01-12
  - 修復 trash 目錄：舊結構移至 `archive/deep/2025-12_old_structure/`
  - 建立二級歸檔結構：`solutions/` 和 `experiments/outputs/` 按任務目標分層
  - 新增 Claude Code Hooks 自動化提醒（程式碼修改、git commit、會話結束）
  - 更新 CLAUDE.md 加入文檔管理規範與 Hooks 說明
  - 建立深度歸檔索引 `archive/deep/2025-12_old_structure/INDEX.md`

- 文檔庫架構重整（第一階段）- 2026-01-11
  - 建立新目錄結構（10 個主要目錄）
  - 遷移 94 個檔案到新架構
  - 統一檔案命名規範（YYYYMMDD_說明_流水號.md）
  - 建立 README.md 文檔庫說明
  - 建立 CURRENT_FOCUS.md 當前目標檔案

## 阻塞問題

（目前無阻塞問題）

---

## 新架構摘要

| 目錄 | 結構 | 用途 |
|------|------|------|
| architecture/ | 扁平 | 專案主軸架構說明 |
| concepts/ | 扁平 | 目標建構過程的構想 |
| plans/ | 年/月 | 迭代計劃表 |
| solutions/ | 任務/年/月 | 問題解決報告 |
| experiments/outputs/ | 任務/年/月 | 測試輸出結果 |
| experiments/parameters/ | 扁平 | 參數研究 |
| ai_sessions/ | 年/月 | AI 對話紀錄 |
| data_specs/ | 類別 | 數據規格 |
| references/ | 類別 | 參考文件 |
| archive/ | 年/月 + deep/ | 歷史歸檔 |

### Hooks 自動化

| 觸發條件 | 提醒內容 |
|----------|----------|
| 編輯 .cpp/.hpp/.h | 編譯 + 測試 |
| git commit | 確認測試 + 文檔 |
| 會話結束 | 撰寫執行報告 |

---

## 快速參考

### 檔案命名格式

```
{YYYYMMDD}_{中文說明目標}_{流水號}.md
```

### 目錄選擇指南

| 我想要... | 放在... |
|-----------|---------|
| 記錄新想法/構想 | `concepts/` |
| 建立正式計劃 | `plans/YYYY/MM/` |
| 記錄問題解決過程 | `solutions/{任務目標}/YYYY/MM/` |
| 保存測試結果 | `experiments/outputs/{任務目標}/YYYY/MM/` |
| 記錄 AI 對話 | `ai_sessions/YYYY/MM/` |
| 查詢歷史結構 | `archive/deep/` |
