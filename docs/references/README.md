<!--
建立時間: 2026-03-12 10:30
目標: 提供 docs/references/ 的總入口，明確拆分 internal manual 與 external references
處理範圍: docs/references/
關聯檔案:
  - docs/README.md
  - docs/references/manual/README.md
  - docs/references/external/README.md
-->

# 參考資料入口

## 目的

`docs/references/` 用來集中管理 **參考性資料**，但必須區分兩種完全不同的用途：

1. `manual/`
   - 專案內部手冊、agent/skill 規格、設定檔與操作說明
2. `external/`
   - 外部文獻、第三方指南、外部方法與整合參考

若不先拆開，使用者容易把「專案內規」和「外部背景」混為同一類資料。

## 子目錄

### 1. 內部手冊與規格

- 入口：
  [docs/references/manual/README.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/README.md)
- 內容：
  AI agent 操作手冊、研究流程手冊、PPT 規範、技能規格、設定檔

### 2. 外部參考

- 入口：
  [docs/references/external/README.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/external/README.md)
- 內容：
  論文整理、第三方工具文件濃縮、外部方法比較

## 邊界規則

1. 專案自己定的規範，不放 `external/`
2. 外部 paper 或第三方教學，不放 `manual/`
3. 若內容已轉化成本專案正式結論：
   - 改放 `docs/reports/` 或 `docs/research/`

## 建議查詢順序

1. 要找專案怎麼做：
   - 先看 `manual/`
2. 要找外部世界怎麼做：
   - 再看 `external/`
3. 要找本專案最後採納什麼：
   - 去 `reports/` 或 `research/`
