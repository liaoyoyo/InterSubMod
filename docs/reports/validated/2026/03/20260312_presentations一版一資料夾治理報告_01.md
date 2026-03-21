<!--
建立時間: 2026-03-12 11:25
目標: 將 docs/presentations/ 收斂為一版一資料夾的輸出治理方式，方便多版本比較、QA 與回查
處理範圍:
  - docs/presentations/validated/2026/03/
  - docs/presentations/draft/2026/03/
  - docs/references/manual/20260311_研究週報PPTX客製化設定與製作手冊_01.md
  - docs/references/manual/assets/20260311_liao_research_ppt_profile_01.json
  - /home/liaoyoyo2001/.codex/skills/liao-research-ppt/SKILL.md
關聯檔案:
  - docs/presentations/README.md
  - docs/reports/validated/2026/03/assets/20260312_presentations_folderization_manifest_01.tsv
-->

# presentations 一版一資料夾治理報告

## 重點結論

`docs/presentations/` 原本雖已有 `draft / validated` 分流，但月目錄下仍直接堆放 `.pptx`、`.version.json`、`pdf/` 與 `assets/`，同一主題多版 deck 會散在不同層級，不利於版本比較、QA 回查與後續再生成。  
本輪已將 `2026-03` 現有 deck 全部整理成 **一版一資料夾**：每個版本資料夾內自帶 `pptx`、`version.json`，需要時再附 `pdf/` 與 `assets/`。同時，`docs/presentations/README.md`、月索引、個人 PPT profile、以及 `liao-research-ppt` skill 都已同步更新，後續新 deck 應直接遵守這個結構。

## 為什麼原本結構不夠好

### 原本的問題

1. 同一月份直接堆多個 `.pptx`
2. `pdf/` 與 `assets/` 放在月目錄共用層
3. 要確認某版 QA 與資產是否完整，必須跨多個位置查
4. `v01 / v02 / v03` 雖然能保留，但不容易一眼看出哪些檔案屬於同一版

### 造成的實際影響

1. 比對版本時容易漏掉對應的 `version.json`
2. 搬動 deck 後，metadata 容易仍指向舊路徑
3. 之後若新增更多版本，月目錄會快速變亂

## 新結構規則

每一版 deck 固定用一個資料夾承接：

```text
docs/presentations/{status}/YYYY/MM/<deck_name>/
├── <deck_name>.pptx
├── <deck_name>.version.json
├── pdf/                    # 可選
└── assets/                 # 可選
```

### 這個結構的好處

1. 同一版的所有產物綁在一起
2. 月目錄只剩 `INDEX.md` 與版本資料夾
3. 可以直接用資料夾作為比較與驗證單位
4. 之後 generator 也更容易直接輸出到固定目標

## 本輪實際調整

### validated

已整理成資料夾：

1. [docs/presentations/validated/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_01/](/big8_disk/liaoyoyo2001/InterSubMod/docs/presentations/validated/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_01/)
2. [docs/presentations/validated/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_oral_v02/](/big8_disk/liaoyoyo2001/InterSubMod/docs/presentations/validated/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_oral_v02/)
3. [docs/presentations/validated/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_oral_v03/](/big8_disk/liaoyoyo2001/InterSubMod/docs/presentations/validated/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_oral_v03/)

### draft

已整理成資料夾：

1. [docs/presentations/draft/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_oral_v01/](/big8_disk/liaoyoyo2001/InterSubMod/docs/presentations/draft/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_oral_v01/)
2. [docs/presentations/draft/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_v04/](/big8_disk/liaoyoyo2001/InterSubMod/docs/presentations/draft/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_v04/)

## 同步更新的規範

1. [docs/presentations/README.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/presentations/README.md)
   - 補上「一版一資料夾」規則
2. [docs/presentations/validated/2026/03/INDEX.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/presentations/validated/2026/03/INDEX.md)
   - 改成列版本資料夾入口
3. [docs/presentations/draft/2026/03/INDEX.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/presentations/draft/2026/03/INDEX.md)
   - 改成列版本資料夾入口
4. [docs/references/manual/assets/20260311_liao_research_ppt_profile_01.json](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/assets/20260311_liao_research_ppt_profile_01.json)
   - 補上 `one_output_per_folder`
5. [docs/references/manual/20260311_研究週報PPTX客製化設定與製作手冊_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260311_研究週報PPTX客製化設定與製作手冊_01.md)
   - 補上輸出資料夾慣例
6. [/home/liaoyoyo2001/.codex/skills/liao-research-ppt/SKILL.md](/home/liaoyoyo2001/.codex/skills/liao-research-ppt/SKILL.md)
   - skill 輸出規則改成一版一資料夾

## 驗證重點

本輪應驗證的不是 deck 內容，而是版本治理是否一致：

1. deck 是否已移入對應資料夾
2. `.version.json` 是否跟 deck 同層
3. `.version.json` 內的 `canonical_deck_path / qa_pdf_path / assets_dir` 是否已指向新位置
4. 月索引是否仍可正確導到 deck
5. 規範檔與 skill 是否已同步寫入新規則

## 變更清單

完整搬移與更新清單見：

[docs/reports/validated/2026/03/assets/20260312_presentations_folderization_manifest_01.tsv](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260312_presentations_folderization_manifest_01.tsv)

## 仍保留但未清掉的項目

由於專案政策不直接刪除檔案，本輪沒有移除月目錄下既有的空 `pdf/` 或空 `assets/` 容器。  
它們目前不再作為主入口使用，但也沒有被刪除。

## 後續建議

1. 後續所有新 deck 都直接輸出到：
   - `docs/presentations/{validated|draft}/YYYY/MM/<deck_name>/`
2. 若 generator 尚未支援此規則，下一輪應直接改生成腳本
3. 之後做版本比較時，以「版本資料夾」為單位，不再以單一 `.pptx` 為單位
