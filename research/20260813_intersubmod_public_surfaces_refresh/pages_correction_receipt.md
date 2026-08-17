<!--
建立時間：2026-08-13
目標：記錄 docs/explain GitHub Pages P1/P2 claim 修正與 SVG/HTML QA。
處理範圍：17 個 *.standalone.html；claim inventory 指定的 24 項 C047–C158 子集。
關聯檔案：InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv
-->

# Pages P1/P2 修正收據

## 結果

- 17/17 個 standalone Pages 已檢查；37/37 個 inline SVG 具有直接子層 `<title>`、`<desc>`、`role="img"` 與有效 `aria-labelledby`。
- page11–16 已補最小 `<!doctype html>`、`<html lang="zh-Hant">`、`<head>`、`<body>` wrapper；未改 SVG 幾何。
- 11 個原有 mobile overflow 的頁面已補 responsive reflow；desktop、390 px mobile、no-JS、print 四種 profile 共 68/68 檢查通過，0 個水平溢出、console error、page error 或 local request failure。
- 8 組公開圖例已由 canonical Pages 來源重新生成 SVG/PNG；8/8 SVG 均為有效 XML 且含非空 `<title>`／`<desc>`，PNG 尺寸介於 810×308 與 978×611。
- `07_subclone-judgment-workstation-chr2-18M.standalone.html` 與其既有 P0 generator 修正均保留，本工作流未重新生成 page07。
- 未修改 README、Wiki、CCU 網站；未 commit、push 或 deploy。

## 補充事實校正

- page11 與 Wiki System Overview 的 HCC1395 tumor BAM 改為 `stat -L` 實測 `283,071,595,503 bytes = 263.63096712436527 GiB`；舊的 `292 GB` 顯示值不再使用。
- 移除沒有可重現 benchmark receipt 的 `2.9 秒`／`約 3 秒`敘述；效能數字須綁定 commit、輸入、硬體、重複次數與分布。
- 移除用 mtime 推斷 binary stale 的說法；現在要求 source commit、rebuild command、toolchain、binary hash 與 test exit code 才能判定 identity。
- host-specific `48 cores / 617 GB` 不再當成目前執行環境；若保留歷史記錄，必須附明確快照日期。
- page15 標題與圖例已限縮為「缺任一已宣告必填 metric 時 fail-closed」，並明示這不會驗證輸入真實性、統計方法或科學結論。

## Claim disposition

| Claim | Pages 處置 |
|---|---|
| C047, C048, C089, C090 | page11/14 改為 exact sidecar `6,256,168,164 bytes (5.826510641724 GiB)`；tagged-BAM 總量、倍率及 `>99%` 欄位占比因缺 receipt/census 不再宣稱。 |
| C110, C111 | page11 將不可識別性限縮為「僅有 per-locus marginal VAF、沒有 linkage/額外假設」及「目前 single-bulk measurement set」。 |
| C114 | index/page11 明列旗標位於 canonical `cohort_receipt.json` 與 `summary/all7_summary.json`，不是 `authority_manifest` top-level。 |
| C115 | 本輪 Pages 無該 README-only claim；不在 Pages 代改。 |
| C116 | 以實際 DOM 元素定義更新/驗證為 37 個 inline SVG。 |
| C117, C118 | index/page01–04/page11 改為「已評估 formulation 的有界負結果」；只有實質新假說加 pre-decision audit 才重啟。 |
| C119 | page01 改為 paired tumor-normal 通常改善區分，但仍須 benchmark、並非精準無誤。 |
| C120, C121 | page01/02 改為 heuristic 下的 cis-compatible / copy-screen-negative，以及 PERMANOVA-positive candidate structure，皆不當作 causal truth。 |
| C122, C123 | page01/02/09/12 明定 p 為與指定 null 的不相容度、V 為關聯幅度；PERMANOVA 僅測 label-associated centroid separation，須併看 PERMDISP 與設計限制。 |
| C132, C158 | index/page02/05 改稱 project-specific combination；novelty/優越性待系統性文獻檢索與 matched benchmark。 |
| C133 | page15 限縮為阻止缺少「已宣告必填指標」時靜默渲染，不保證科學正確性。 |
| C145 | index/page11/16 更新為 HEAD `95d420f6f397`、270 tests / 39 suites、7039 ms、exit 0。 |
| C146, C147 | page15 更新 deploy commit `fbdf7c7` 為全 repo 2,147 個 `.py`、`scripts/` 291 個 `.py`，並刊出計數命令及包含範圍。 |
| C156 | page10 改為 matched normal 沒有六個預期 ALT，但有零星其他 non-reference errors。 |
| C157 | Pages 目前沒有 `32 cores linear speedup / <300 ms` claim；不在 Pages 補造 benchmark。 |

## 驗證收據

1. 輸入：`InterSubMod/docs/explain/*.standalone.html`  
   命令：`python3 tools/explain_page_qa.py docs/explain/*.standalone.html`  
   輸出：stdout；17 頁逐頁通過，`總計 0 個 FAIL`（無 SVG 的 3 頁僅為資訊性 WARN）。
2. 輸入：同上 17 頁  
   命令：BeautifulSoup 結構、ID、SVG title/desc/ARIA 與內部 HTML link audit  
   輸出：`pages=17 svg=37 svg_pass=37 svg_fail=0 issues=0`，exit 0。
3. 輸入：`InterSubMod/build/bin/run_tests`  
   命令：`./build/bin/run_tests`  
   輸出：270 tests / 39 suites passed，7039 ms，exit 0；HEAD `95d420f6f397`。
4. 輸入：全部 Pages 文字  
   命令：`rg` guard 搜尋舊的 1.67 TiB、287×、>99%、265/38、2,144、268、32 cores/<300 ms 與過度因果措辭  
   輸出：0 筆，exit 0（以 `|| true` 呈現空結果）。
5. 輸入：`InterSubMod/docs/explain/*.standalone.html`  
   命令：本地 HTTP server + `python3 research/20260813_intersubmod_public_surfaces_refresh/scripts/browser_qa_explain_pages.py --base http://127.0.0.1:8765/docs/explain --output research/20260813_intersubmod_public_surfaces_refresh/pages_browser_qa.json`  
   輸出：17 pages × 4 profiles = 68 checks；`failures=[]`、`verdict=PASS`，所有頁面 HTTP 200 且無水平 overflow。
6. 輸入：`InterSubMod/tools/extract_svg_for_github.py` 與 `InterSubMod/docs/explain/*.standalone.html`  
   命令：`python3 tools/extract_svg_for_github.py`，再以 `xmllint --noout docs/images/*.svg` 與 PIL 讀取全部 PNG  
   輸出：8 組 SVG/PNG；8/8 XML PASS、8/8 title/desc PASS、8/8 PNG 可讀；舊 `2.9 秒`、`292 GB` 與過度 fail-closed claim guard 0 hit。
