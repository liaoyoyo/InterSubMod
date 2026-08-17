<!--
建立時間：2026-08-13
目標：保存 InterSubMod 公開資訊與 CCU 清單技術報告的 canonical artifact、HTML packaging、SVG extraction 與 browser QA 證據
處理範圍：本輪唯一 primary report surface；不含 GitHub/Pages/CCU publication
關聯檔案：
  - InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/report_artifact.json
  - InterSubMod/docs/reports/validated/2026/08/20260813_InterSubMod公開資訊更新與CCU改進清單_01.standalone.html
-->

# HTML report build receipt

## 結論

Canonical report artifact 與 self-contained HTML 已完成。最終 HTML 由 Data Analytics portable
report 的官方 validator、reader、chart SVG extractor 與 browser verifier 建置／驗證；包含 25 個
manifest blocks，reader 展開為 30 個 rendered content blocks、3 張 native charts（light/dark 共
6 個內嵌 static SVG）、4 張 tables 與 6 張 metric cards。1440 px desktop 與
390 px narrow viewport 均通過，source dialog 與 keyboard interaction 通過。

唯一 renderer compatibility 例外是：官方 shared reader 在 headless-shell iframe 以 `100vw` 計算
頂欄，而 classic scrollbar 會使 `innerWidth` 比 `clientWidth` 多 15 px。最終 adapter 只加入
`scrollbar-gutter:stable` 與 content-width 頂欄規則，再呼叫相同官方 verifier；沒有改 canonical
artifact payload、沒有加入第二套 chart runtime，也沒有用 `overflow:hidden` 隱藏錯誤。

## 輸入與輸出

### 輸入

- `InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/report_artifact.json`
- `InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/scripts/deliver_public_refresh_report.mjs`
- Data Analytics plugin root：版本路徑 `data-analytics/0.2.8-13ceeea1f599`

### 輸出

- Primary HTML：`InterSubMod/docs/reports/validated/2026/08/20260813_InterSubMod公開資訊更新與CCU改進清單_01.standalone.html`
- Supporting Markdown：`InterSubMod/docs/reports/validated/2026/08/20260813_InterSubMod公開資訊更新與CCU改進清單_01.md`
- Canonical artifact：`InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/report_artifact.json`
- Failure-only diagnostic：`InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/report_packaging_failure_desktop.png`

## Step → Verify

1. 由 inventory／registries／receipts 產生 artifact → 驗證：腳本 assert `158 claims`、`58 problems`、
   `34/20/4 priorities` 與兩個 registry exact-set closure。
2. 驗證 provenance contract → 驗證：12 個 snapshot datasets 各有實際
   `SELECT * FROM <dataset>;` query source；全部檔案來源為 repo-relative path，0 absolute／parent traversal。
3. 以官方 portable renderer 建 HTML 與 static SVG → 驗證：artifact validation／package／
   3-chart SVG extraction 全 PASS。
4. 套用 scrollbar compatibility rule → 驗證：官方 browser verifier 在 1440/390 兩 viewport
   均無 horizontal overflow，counts、geometry、source dialog、keyboard interaction、external request、
   console/page error 全 PASS。
5. 保存 exact hash → 驗證：MD、HTML、artifact、CCU checklist 均有 SHA-256，供 publication 後比對。

## 執行命令與輸出

### Artifact builder

```bash
python3 research/20260813_intersubmod_public_surfaces_refresh/scripts/build_public_refresh_report_artifact.py
```

實際輸出：

```text
artifact=research/20260813_intersubmod_public_surfaces_refresh/report_artifact.json
claims=158 problem=58 dispositions=58 pages=17 svg=37 browser_checks=68 ccu_remaining=16
```

### Final portable delivery

```bash
node research/20260813_intersubmod_public_surfaces_refresh/scripts/deliver_public_refresh_report.mjs \
  --plugin-root /bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599 \
  --input research/20260813_intersubmod_public_surfaces_refresh/report_artifact.json \
  --output docs/reports/validated/2026/08/20260813_InterSubMod公開資訊更新與CCU改進清單_01.standalone.html
```

實際輸出摘要：

```json
{
  "ok": true,
  "officialStages": {
    "validationAndPackage": "passed",
    "staticChartExtraction": "passed",
    "verification": "passed"
  },
  "verification": {
    "counts": {"blocks": 30, "charts": 3, "html": 0, "metrics": 6, "tables": 4},
    "sourceDialog": "passed",
    "sourceInteraction": "keyboard_menu_semantic_click",
    "viewports": [1440, 390]
  }
}
```

### Final hashes

```text
99b4fc5708d8c44eec46eb9d7cc41f3918679bb666288ca582b9bf7d96c8c411  supporting Markdown
41e615706e6030eb6f76aa4d07ce823858ef9d22f9b8a45d53d3e8139e7a2b65  primary HTML
b03365273d7e58bed4bbed73836c12b1452bd7df0ef158044bbb631d8d8d593c  canonical artifact
b19eae989d4c3db8748a5f9e93f162d8396d5ac7a6b8e8f282e71d5b42d14dd2  CCU checklist
```

## 失敗收據與修正

1. 第一次官方 delivery 在 package gate 拒絕：cards/charts/tables 的 source 缺實際 SQL text。
   修正為每個 snapshot dataset 增加可讀的 `SELECT * FROM <dataset>;` query metadata，沒有放寬 validator。
2. 第二次官方 delivery 在 desktop iframe 回報 horizontal overflow；diagnostic 證明表格都在
   `overflow-x:auto` 容器內，未裁切元素清單為空，問題是 top bar `100vw`＋classic scrollbar。
3. 加入 content-width 規則後，mobile iframe 揭露相同 scrollbar geometry；增加
   `scrollbar-gutter:stable` 後，官方 verifier 兩個 viewport 全部 PASS。

第一次 desktop failure-only screenshot 已移到 research cycle 保留，沒有刪除或冒充成功輸出。

## 最終判定

**PASS — self-contained HTML report delivered and browser-verified.** 這個判定只涵蓋報告 artifact、
SVG、tables、responsive reader 與 source affordances；不表示 GitHub/Pages 已發布，也不表示 CCU 已修復。
