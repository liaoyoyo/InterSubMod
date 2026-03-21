<!--
建立時間: 2026-03-11 20:10
更新時間: 2026-03-12 11:20
目標: 提供 InterSubMod 簡報輸出的單獨入口、版本確認方式、draft/validated 分流與一版一資料夾查找規則
-->

# InterSubMod 簡報入口

## 目的

`docs/presentations/` 是 `.pptx` / `.pdf` 簡報的專用入口。  
研究報告的 Markdown 維持放在 `docs/reports/`，但簡報另外有單獨資料夾，方便：

1. 直接檢視 deck
2. 確認版本與 QA 狀態
3. 對照來源報告
4. 查找生成腳本、設定檔與 deck 資產

## 目錄分流規則

1. `validated/`
   - `docs/presentations/validated/YYYY/MM/`
   - 只放已完成 QA、可對內引用或可交付的 deck
2. `draft/`
   - `docs/presentations/draft/YYYY/MM/`
   - 放尚未完成 QA 的 deck、較早期草稿、或已被新版取代但仍需保留的設計迭代
3. 月目錄內採 **一版一資料夾**
   - `docs/presentations/{status}/YYYY/MM/<deck_name>/`
   - 月目錄本身只留 `INDEX.md` 與版本資料夾

## 版本治理最低要求

### validated deck

每份 `validated` 輸出資料夾至少必須包含：

1. deck 檔（`.pptx`）
2. 版本確認檔（`.version.json`）
3. 來源報告路徑
4. 生成腳本路徑
5. profile 或設定檔路徑
6. QA 狀態
7. 若為正式交付版，建議同時保留 `pdf/`
8. 若有圖資或渲染輸出，集中放在 `assets/` 或其他子目錄

### draft deck

每份 `draft` 輸出資料夾至少應包含：

1. deck 檔（`.pptx`）
2. 版本確認檔（`.version.json`）
3. 狀態標記（例如 `draft`、`superseded_draft`、`qa_pending`）
4. 來源報告路徑
5. 若有資產或設定檔，應在 `.version.json` 中指向對應路徑

## 輸出資料夾格式

```text
docs/presentations/{status}/YYYY/MM/<deck_name>/
├── <deck_name>.pptx
├── <deck_name>.version.json
├── pdf/                    # 可選，正式傳閱版
└── assets/                 # 可選，圖資、渲染、插圖、QA 輔助資產
```

這樣做的目的：

1. 同一版的 `pptx / version / pdf / assets` 綁在一起
2. 同一月份可並排保留多個版本資料夾，方便比較與 QA
3. 月目錄不再直接堆滿 deck 檔

## deck 家族與命名

1. `*_oral_vXX.pptx`：口頭報告型 deck
2. `*_vXX.pptx`：設計迭代 deck 或特定版型 deck
3. `*.version.json`：版本、QA、來源與重生設定
4. 每個版本各自一個資料夾，不共用同層 `pdf/` 與 `assets/`
5. `pdf/`：交付或傳閱用 PDF
6. `assets/`：圖表、插圖與生成資產

## 目前入口

1. [2026-03 validated PPT index](/big8_disk/liaoyoyo2001/InterSubMod/docs/presentations/validated/2026/03/INDEX.md)
2. [2026-03 draft PPT index](/big8_disk/liaoyoyo2001/InterSubMod/docs/presentations/draft/2026/03/INDEX.md)
