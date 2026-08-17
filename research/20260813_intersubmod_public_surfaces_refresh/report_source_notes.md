<!--
建立時間：2026-08-13
目標：保存公開資訊校正技術報告的 job definition、section mapping、chart map、來源與省略理由
處理範圍：HTML technical report；不修改或發布 GitHub／CCU external surfaces
關聯檔案：
  - InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/report_artifact.json
  - InterSubMod/docs/reports/validated/2026/08/20260813_InterSubMod公開資訊更新與CCU改進清單_01.standalone.html
-->

# Report source notes — InterSubMod public surfaces refresh

## Reporting job

- **Question**：InterSubMod GitHub／Pages 哪些說明需要修正、目前修到哪裡、哪些結論可保留；CCU 還需要補什麼？
- **Audience**：technical，重點是 claim/data/method/publication validation。
- **Scope**：158 個 frozen claim families、58 個 problem claims、17 個 Pages、37 個 inline SVG、四種 browser profiles；CCU 32 個 frozen findings。
- **Time frame**：evidence frozen 2026-08-12 至 2026-08-13；live GitHub/Pages state 於 2026-08-13 查詢。
- **Comparison basis**：公開敘述對 authority artifacts、code/runtime receipts、claim registries 與 live refs/bytes；Pages accessibility 比較修正前後同一組 37 SVG。
- **Decision-useful success**：使用者可以分辨 local corrected、live resolved、publication pending 與 checklist-only；可以看到明確的修正清單、claim ceiling 與下一步 gate。

## Technical report required-structure mapping

| Technical specification role | Visible report section |
|---|---|
| Title | `InterSubMod 公開資訊校正與 CCU 改進稽核` |
| Technical summary | `技術摘要 — 本地公開來源已校正，live 發布仍未完成` |
| Key findings with visual evidence | claim verdict、SVG accessibility、CCU partition 三段 |
| Scope, data, metric definitions | `推論上限已統一為 molecule → local candidate → association`，各 chart denominator subtitle |
| Methodology | `驗證方法 — inventory、registry、負控制、結構與 live bytes 分層` |
| Limitations / uncertainty / robustness | `限制與 robustness — 文件正確不等於方法已被 truth-set 驗證` |
| Recommended next steps | `建議下一步 — 先發布 InterSubMod，再另行決定是否授權 CCU 修復` |
| Further questions | `仍會改變決策的問題` |

沒有省略 required role。Scope/definitions 在 key findings 後完整呈現，但每張圖的 grain、denominator 與 unit 已在圖表前／subtitle 定義，避免讀者先看到無分母數字。

## Chart map

| Section | Analytical question | Family/type | Fields | Supported claim | Palette / non-color | Final surface |
|---|---|---|---|---|---|---|
| Claim audit | 158 個 claim 各 verdict 有多少？ | Comparison / bar | verdict → claims | 58 個 problem claims 是 28+26+4；31 項 limits 不應誤算漏修 | 單系列；axis labels + exact tooltip，不以顏色製造群組 | portable HTML native SVG |
| Pages QA | 37 個 SVG 的 accessibility contract 是否完整閉合？ | Group comparison / bar | stage × status → svg | 修正前 11 pass/26 missing；修正後 37 pass/0 missing | status 用兩類色；同時保留 stage/status 文字與數字 | portable HTML native SVG |
| CCU closure | 32 個 finding 的三種 disposition 各多少？ | Comparison / bar | partition → findings | 13 patch-only、16 remaining、3 prior-resolved；patch 不等於 live repair | 單系列；類別文字與 exact tooltip | portable HTML native SVG |

重複使用 bar 的理由：三段問題都是少量離散類別的 exact count comparison，沒有時間序列、分布或連續關係；改用 line、pie、scatter 會降低可讀性或暗示不存在的關係。

## Evidence and calculation notes

- Claim verdict 與 priority 直接從 `claim_inventory.tsv` 重算；builder assert `158`、`69/31/28/26/4`、problem `58`、priority `34/20/4`。
- Disposition closure 以 P0 與 P1/P2 registry claim-ID 聯集對 problem set 做 exact-set equality；不是以文字命中數估計。
- SVG 修正前 26 missing 來自本 cycle implementation notes 的 frozen initial QA；修正後 37/37 來自 canonical structural QA JSON。
- Browser `68=17 pages×4 profiles`；final receipt `failures=[]`。早期 harness/overflow failures 留作診斷，不混入 final pass dataset。
- CCU `32=13+16+3` 由 frozen TSV ID partition 驗算；16 remaining 的七個 acceptance 欄位均 16/16。
- Public live status只採 remote receipt；本地 git diff 不用來宣稱遠端已更新。

## Omitted visuals and metrics

- 不畫「1.67 TiB → 5.83 GiB」縮減圖：tagged-BAM total 缺 exact path/byte/hash/compression receipt，兩端不是可比對象。
- 不畫效能/speedup 圖：`2.9 秒`、`<300 ms`、`32-core linear speedup` 缺 commit/input/hardware/repetition/distribution receipt，公開數值已撤回。
- 不畫 F1 或 truth-set improvement：本輪是 documentation/publication audit，沒有新增 caller-matched TP/FP/FN benchmark。
- 不把 3/811 甲基化 association 畫成方法效能：這是 frozen yield，不是 sensitivity、specificity 或所有可能甲基訊號的總檢定。

## Report delivery boundary

唯一 primary surface 是自包含 `.standalone.html`；`report_artifact.json` 是 canonical renderer input，validated Markdown 是 supporting technical source。未建立平行 dashboard、MCP report、Sites deployment 或 PDF。

### Renderer compatibility note

官方 `report:deliver` 已通過 artifact/package 與 SVG extraction，但在 headless-shell 的 desktop iframe 因 shared `.analytics-top-bar { width:100vw }` 加上 classic vertical scrollbar，多出 15 px document overflow。針對性診斷證明長表位於自己的 `.table-wrap { overflow-x:auto }`，沒有未裁切元素越過 viewport；問題是 top bar 的 viewport-width 計算。

本輪因此使用 `deliver_public_refresh_report.mjs` 呼叫同一 plugin 的官方 `buildPortableArtifact`、`extractPortableChartSvgs` 與 `verifyPortableArtifact`，只在 self-contained HTML 加入 `scrollbar-gutter:stable` 與 scoped top-bar/fallback-header content-width rule，再由官方 verifier 重驗。沒有加入第二套 chart/runtime，沒有用 `overflow:hidden` 隱藏失敗，也沒有更改 canonical artifact payload。
