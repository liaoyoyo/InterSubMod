<!--
建立時間: 2026-08-13 Asia/Taipei
目標: 保存 P0 修正 standalone HTML 的桌機、手機、列印、連結、SVG 與視覺 QA 證據
處理範圍: InterSubMod/docs/reports/validated/2026/08/20260813_公開文件P0修正與驗證_01.standalone.html
關聯檔案:
  - InterSubMod/research/20260813_public_docs_p0_correction/report_qa/qa_receipt.json
  - InterSubMod/research/20260813_public_docs_p0_correction/scripts/qa_correction_report.py
-->

# Standalone HTML／SVG 終局 QA 收據

用 PREP：**fresh-reader 修訂後的 browser QA 通過 — 桌機／手機 0 overflow，0 console/page error，0 external request，4/4 SVG 可存取，列印時 7 個折疊內容全可見，exit 0**（影響：中；信心：高）。

## 輸入、命令與輸出

- 輸入：`InterSubMod/docs/reports/validated/2026/08/20260813_公開文件P0修正與驗證_01.standalone.html`
- 命令：

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 research/20260813_public_docs_p0_correction/scripts/qa_correction_report.py \
  docs/reports/validated/2026/08/20260813_公開文件P0修正與驗證_01.standalone.html \
  --out-dir research/20260813_public_docs_p0_correction/report_qa
```

- 輸出：
  - `InterSubMod/research/20260813_public_docs_p0_correction/report_qa/qa_receipt.json`
  - `InterSubMod/research/20260813_public_docs_p0_correction/report_qa/desktop.png`
  - `InterSubMod/research/20260813_public_docs_p0_correction/report_qa/mobile.png`
  - `InterSubMod/research/20260813_public_docs_p0_correction/report_qa/print.pdf`

## 迭代執行與修正

第一次執行的內容 gate 均通過，但 QA 將 `left:-9999px` 的 skip-link 視為 horizontal overflow，故 `pass=false`、exit 1。這是 accessibility CSS pattern 與 checker 的衝突，不是內容或 SVG failure。

修正：改成 width/height 1px + `clip-path: inset(50%)` 的 visually-hidden pattern，focus 時恢復尺寸；未移除 skip-link。

獨立 fresh-reader 後補入 uncommitted manifest 說明時，完整 SHA-256 與長 badge 令手機 `documentScroll=535`，QA 再次以 exit 1 阻擋。修正為可換行 footer 與短識別，完整 hash 保留在 manifest；重跑後 `390/390`。同輪也修正圖 1 五段 bar 的幾何比例，消除原先重疊。

## 最終實際輸出

```text
size_bytes=40410
desktop: viewport=1440, documentScroll=1440, uncontainedWide=[]
mobile: viewport=390, documentScroll=390, uncontainedWide=[]
h1=1, h2=9, tocLinks=9, details=7, tables=4
console_errors=[]; page_errors=[]; external_requests=[]
broken_internal_links=[]; missing_anchor_targets=[]
svg=4; all viewBox/role/title/desc/aria-labelledby=true
hidden_print_detail_bodies=0
gradient_count=0; external_asset_markup=false; oversize=false
pass=true
EXIT=0
```

## 人工視覺檢查

- 桌機：sticky TOC、九段敘述、表格、callout 與四張 SVG 層級清楚；圖 1 五段互不重疊，沒有長 SHA 逸出。
- 手機：單欄閱讀、stats 疊放、表格在其容器內橫向捲動；整頁無水平裁切。
- 列印：A4 PDF 成功，closed details 的 `.card-body` 仍顯示。

![桌機全頁 QA](report_qa/desktop.png)

![手機全頁 QA](report_qa/mobile.png)

## 判定邊界

本收據證明報告 artifact 的呈現、連結、runtime dependency 與 SVG accessibility gate；不證明底層 biological claim，也不代表 GitHub／CCU live site 已部署。
