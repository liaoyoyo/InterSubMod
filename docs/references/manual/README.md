<!--
建立時間: 2026-03-12 10:30
目標: 提供 docs/references/manual/ 的手冊分類、查詢順序與 assets 邊界
處理範圍: docs/references/manual/
關聯檔案:
  - docs/references/README.md
  - docs/CURRENT_FOCUS.md
  - docs/experiments/INDEX.md
-->

# 內部手冊與規格入口

## 目的

`docs/references/manual/` 放的是 **專案內部操作手冊、agent/skill 規格與設定檔**。  
它不是正式研究結論入口，而是告訴你：

1. 每次研究前要讀什麼
2. 如何記錄、驗證與整理研究
3. AI agent / skills 如何使用
4. PPT 與文件如何依個人風格輸出

## 建議查詢順序

### A. 研究啟動

1. [docs/CURRENT_FOCUS.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md)
2. [docs/experiments/INDEX.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/INDEX.md)
3. [docs/references/manual/20260307_研究推進與實驗觀察手冊_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260307_研究推進與實驗觀察手冊_01.md)

### B. 文件 / Agent / Skills

1. [docs/references/manual/20260310_InterSubMod研究文件Agent規格_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_InterSubMod研究文件Agent規格_01.md)
2. [docs/references/manual/20260310_研究報告Agent與Skills使用手冊_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_研究報告Agent與Skills使用手冊_01.md)
3. [docs/references/manual/20260310_研究脈絡整理Skill規格_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_研究脈絡整理Skill規格_01.md)
4. [docs/references/manual/20260310_週報生成Skill規格_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_週報生成Skill規格_01.md)
5. [docs/references/manual/20260310_指令修正與偏好收斂Skill規格_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_指令修正與偏好收斂Skill規格_01.md)

### C. PPT 與個人風格

1. [docs/references/manual/20260311_個人PPT設計風格規範_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260311_個人PPT設計風格規範_01.md)
2. [docs/references/manual/20260311_研究週報PPTX客製化設定與製作手冊_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260311_研究週報PPTX客製化設定與製作手冊_01.md)
3. `assets/` 下的 JSON profile 與 config

## 子類別

1. 研究流程與操作手冊
2. Agent / Skill 規格
3. 個人風格規範與設定檔
4. 小型歷史 checklist / 速查筆記

## assets 邊界

`docs/references/manual/assets/` 只放：

1. JSON 設定檔
2. profile
3. 可機器讀取的生成 config

這裡不放：

1. 簡報圖資
2. 研究結果圖
3. 大型輸出

## 現有 assets

1. [docs/references/manual/assets/20260311_liao_research_ppt_profile_01.json](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/assets/20260311_liao_research_ppt_profile_01.json)
2. [docs/references/manual/assets/20260311_weekly_report_pptx_config_01.json](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/assets/20260311_weekly_report_pptx_config_01.json)
3. [docs/references/manual/assets/20260311_weekly_report_pptx_oral_config_01.json](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/assets/20260311_weekly_report_pptx_oral_config_01.json)
4. [docs/references/manual/assets/20260311_weekly_report_pptx_oral_config_02.json](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/assets/20260311_weekly_report_pptx_oral_config_02.json)
5. [docs/references/manual/assets/20260311_weekly_report_pptx_oral_config_03.json](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/assets/20260311_weekly_report_pptx_oral_config_03.json)
