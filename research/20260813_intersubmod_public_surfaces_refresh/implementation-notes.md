<!--
建立時間: 2026-08-13 17:36 +08:00
目標: 即時記錄 InterSubMod 公開介面資訊更新的設計決定、偏離、折衷與未決事項
處理範圍: 本輪 GitHub/Pages 修改、CCU 唯讀清單與驗證鏈
關聯檔案:
  - InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/00_INDEX.md
  - InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/pre-decision-audit.md
-->

# Implementation notes — InterSubMod public surfaces refresh

狀態：`validated_local / publication_pending`  
開始時間：2026-08-13 17:36 +08:00  
Task type：B — Comprehensive validation

## 設計決定

### D-001 — 使用問題 verdict 與 priority 的交集凍結 24 claims

- 決定：只有 verdict=`NEEDS_CORRECTION|CONTRADICTED|UNVERIFIABLE` 且 priority=`P1|P2` 才列入本輪 correction denominator。
- 理由：inventory 內另有 31 個 `CONFIRMED_WITH_LIMITS/P2`；它們是保留限制的已確認項，不能誤算成待修問題。
- 驗證：六類計數為 `P1: 6 contradicted + 11 needs correction + 3 unverifiable`、`P2: 3 needs correction + 1 unverifiable`，合計 24。

### D-002 — scientific claim ceiling

- 直接觀測上限：同一條 physical molecule 上的 allele/HP/PS/MM/ML called evidence。
- 結構上限：local recurrence-allowed model-conditional mutation-state candidates。
- 甲基上限：genetic-pattern-conditioned regional methylation association。
- 禁止升格：confirmed cellular clone/lineage、canonical CN/LOH-corrected CCF、causal methylation/function。

### D-003 — 使用者保護決策：CCU 僅清單

- **[決策／不可偏離] 使用者明示「CCU 教學站只整理重點改進清單，不改動 CCU」。**
- 允許：唯讀既有 audit/receipt，在本 repo 新增 CCU 改進清單。
- 禁止：改 CCU source、remote、live site、既有 CCU patch，或產生新的可套用 CCU patch。
- fail condition：任何 CCU material checksum 或 git state 因本輪改變。

### D-004 — 無 receipt 的精確數字採 fail-closed

- tagged-BAM `1.67 TiB`、`287×`、SEQ/QUAL `>99%` 與 Project Summary performance 數字，在缺 exact receipt／field census／benchmark 時不以估計值替換，而是移除或標 `UNVERIFIED`。
- 可保留的 exact 值：current 7 sidecars 合計 `6,256,168,164 bytes = 5.826510641724 GiB`，其重算命令已在 2026-08-12 command receipt 保存。

### D-005 — 版本統計必須綁 corpus

- `2,147` tracked `.py`、`scripts/` 內 `291` 個 tracked `.py`、37 個 inline SVG elements 都只適用於 deploy `fbdf7c7`，並附 `git ls-tree`／`git grep` counting command。
- `270 tests / 39 suites` 只屬 tracked core `73afaeac`、2026-08-12 snapshot；未來版本以 fresh command output/exit 為準。

### D-006 — 公開頁面需同時通過內容與呈現驗證

- 內容 gate：P0 與 P1/P2 claim guard 都必須 PASS，且負控制須能偵測遺漏 claim。
- 結構 gate：17/17 HTML、37/37 inline SVG title/desc/ARIA 與 local links 必須 PASS。
- 呈現 gate：desktop、390 px mobile、no-JS、print 共 68/68 browser checks 必須無 overflow、console/page error 或 local request failure。
- 可攜性 gate：QA JSON 與最終 artifact 不得暴露工作站絕對路徑。

## 偏離

### V-001 — C108 GitHub About 已在 live state 解決

- 2026-08-13 17:49 +08:00 重新查詢 GitHub API 與 repository page，About 已是有界新描述。
- 這是 2026-08-13 P0 receipt 之後的外部狀態變更；新報告將 C108 更新為 `RESOLVED_LIVE`，不回寫或改動 immutable 舊 audit。
- Default `main` 與 17 個 Pages 尚未變更，因此其他 remote publication status 不隨 C108 一起升格。

### V-002 — GitHub public-source claim set 出現一個 derived Wiki drift

- Frozen inventory 的 GitHub occurrence 為 14 個問題 claims；current `docs/wiki/System-Overview.md` 另有 C117 `DEAD／勿再開` 的未登記 recurrence。
- 處置：在不變更既定 24-claim denominator 下修正該 derived drift；GitHub source 結果為 14 個既定 occurrence claims + 1 個 current drift，剩餘 9 個 claims 為 Pages-only。

### V-003 — CCU 清單終局驗證

- `32 = 13 prior patch targets + 16 remaining + 3 prior resolved`，unclassified=0、duplicates=0。
- 16 個 remaining 的七個必要欄位均為 16/16；既有 patch SHA-256 仍為 `aee944790ad8b60634818fca8888b0f980b3c1e5522ef64c00015a0a21d6c3e6`。
- 本輪 CCU disposition 固定為 `CHECKLIST_ONLY / NO_CCU_CHANGE / NOT_PROOF_OF_REPAIR`。

### V-004 — Pages 響應式與圖例終局驗證

- browser QA 曾依序揭露 harness base mismatch、14、5、2 個 overflow failures；各輪失敗收據保留作診斷，但只把最終 `pages_browser_qa.json` 列為 canonical PASS receipt。
- 最終結果：17 pages × desktop/mobile/no-JS/print = 68/68 checks PASS；`failures=[]`。
- 37/37 inline SVG 結構 PASS；8 組 GitHub 圖例重新生成後，SVG XML/title/desc 與 PNG decoding 8/8 PASS。
- HCC1395 tumor BAM 以 `stat -L` 實測為 `283,071,595,503 bytes = 263.63096712436527 GiB`；舊 `292 GB` 已撤回。

### V-005 — 主交付 HTML 封版

- canonical artifact 含 25 個 manifest blocks；official reader 實際呈現 30 blocks、3 charts、4 tables、6 metrics。
- 三幅 chart 均有 native SVG 靜態圖；日間／夜間圖合計嵌入 6 個 `<svg>`。
- official portable-artifact 三階段均 PASS：`validationAndPackage`、`staticChartExtraction`、`verification`；1440 px 與 390 px 的 source dialog 與 keyboard semantic click 均 PASS。
- 封版 SHA-256：Markdown `99b4fc5708d8c44eec46eb9d7cc41f3918679bb666288ca582b9bf7d96c8c411`；HTML `41e615706e6030eb6f76aa4d07ce823858ef9d22f9b8a45d53d3e8139e7a2b65`；artifact `b03365273d7e58bed4bbed73836c12b1452bd7df0ef158044bbb631d8d8d593c`。

### V-006 — 最終回歸關閉

- P0 claim guard：34/34 PASS，`errors=0`。
- P1/P2 claim guard：24/24 PASS，`errors=0`；C157 負控制正確 exit 1。
- Pages：17/17 頁與 37/37 inline SVG 結構 PASS；68/68 browser checks PASS。
- Core regression：270 tests / 39 suites PASS。
- 上述只證明文件、呈現與回歸契約無已知缺口，不等於新的生物真實性或效能驗證。

### V-007 — evidence ledger 可追蹤性修復

- 本輪只 append 第 151 筆 entry，未改寫前 150 筆。
- `research/AGENTS.md` 將 `research/autoresearch/evidence_ledger.jsonl` 定義為 append-only SoT，但廣域 `research/**/*.jsonl` 使它被 Git 忽略；本輪在 `.gitignore` 增加單檔精確 negation，不改其他 JSONL 政策。
- JSONL parse PASS：151/151，last cycle 為 `20260813_intersubmod_public_surfaces_refresh`。

## 折衷

### T-001 — 本地來源完整修正，但 remote publication 分離

- 本輪可完整修正與驗證 local source；GitHub About 已在 live state 解決，default `main`、Wiki 與 Pages live bytes 仍需另一次具發布權限的操作。
- 報告必須分開標示 `LOCAL_SOURCE_CORRECTED` 與 `LIVE_PUBLICATION_PENDING`，不能把未發布狀態說成已上線。

## 未決

- U-001：7 個 tagged-BAM 的 exact paths/bytes/compression/hash receipt 未取得；相關 total/ratio 保留 `UNVERIFIED`，不得用 raw BAM 或 rounded display 值替代。
- U-002：C157 performance claim 缺 hardware/commit/input/repetitions/distribution/scaling artifact；在可重現 benchmark 出現前移除公開數字。
- U-003：GitHub About 已由 live API 證實修正；default `main`、Wiki 與 Pages 仍是舊 bytes，必須在另一次明確授權的 publication 動作後再驗證。

## Lore

- 2026-08-13 P0 correction 已完成本機 33 個 claim，C108 GitHub About 為 external action；P0 guard 34/34 PASS。
- 既有 Pages QA：17 頁、37 SVG、0 structural FAIL；26 SVG 缺 `title`/`desc`，分布於 02/03/04/05/06/08/09。
- CCU 已知 finding closure：32 = prior patch targets 13 + remaining unresolved 16 + prior resolved 3；本輪只處理資訊整理。
