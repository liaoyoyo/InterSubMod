<!--
建立時間: 2026-03-12 00:40
目標: 定義 docs/references/external/ 的收錄範圍、命名方式與引用規則
-->

# 外部參考入口

## 目的

`docs/references/external/` 專門收錄：

1. 外部文獻整理
2. 第三方工具指南整理
3. 外部方法、標準或教學的濃縮筆記
4. 與專案整合相關的外部方案整理

這裡**不放**：

1. 專案內部操作手冊
2. AI workflow、agent 與 skill 規格
3. 專案自己的正式研究結論

## 結構規則

1. 年月分層：
   - `docs/references/external/YYYY/MM/`
2. 檔名：
   - `YYYYMMDD_主題_流水號.md`
3. 每份文件都應補 metadata，並標出：
   - 外部來源主題
   - 與 InterSubMod 的關聯
   - 是否已驗證或只作背景參考

## 建議內容欄位

1. 外部來源類型：
   - paper
   - vendor docs
   - method review
   - integration guide
2. 為何需要這份參考
3. 與本專案直接相關的結論
4. 不確定或待驗證之處
5. 原始來源連結或可追溯資訊

## 與其他 docs 類別的邊界

1. `docs/references/manual/`
   - 專案內部手冊、agent/skill 規格、設定檔
2. `docs/research/`
   - 本專案自己的研究脈絡與研究文檔
3. `docs/reports/`
   - 本專案正式驗證報告與結論

## 現有文件

1. [docs/references/external/2026/01/20260122_Antigravity_ClaudeCode_Integration_Guide_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/external/2026/01/20260122_Antigravity_ClaudeCode_Integration_Guide_01.md)
2. [docs/references/external/2026/03/20260302_tissue_origin_methylation_confounding_literature_review_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/external/2026/03/20260302_tissue_origin_methylation_confounding_literature_review_01.md)
